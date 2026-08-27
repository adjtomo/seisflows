#!/usr/bin/env python3
"""
The Portable Batch System (PBS), and its derivatives (PBS Professional,
OpenPBS, TORQUE) is a widely used workload manager on HPC systems. The Pbs
system class provides generalized utilities for interacting with PBS-based
systems.

.. note::
    PBS deployments differ quite a bit between sites/vendors (OpenPBS, PBS
    Professional, and further site-specific customization on top of either).
    This module targets behavior that is common to PBS-derivatives, i.e.,
    job submission through 'qsub' using '#PBS' directives, job monitoring
    through 'qstat', and job deletion through 'qdel'. Job IDs are assumed to
    take the PBS Professional form '<sequence number>[.<index>].<server>'.

    Because job STATUS values reported by 'qstat' (e.g., 'R', 'Q', 'F')
    typically do not by themselves distinguish a successfully completed job
    from a failed one, this module additionally relies on 'tracejob' to
    inspect a finished job's 'Exit_status' attribute (0 == success). If a
    given PBS deployment does not provide 'tracejob', or reports states in a
    different vocabulary, `query_job_states` will need to be overwritten by
    the relevant child class (see e.g., `system.Miyabi`).

.. note::
    Unlike SLURM's 'sbatch', PBS's 'qsub' does not reliably forward
    trailing command line arguments to the submitted script -- some 'qsub'
    builds scan the *entire* command line for dash-prefixed tokens (GNU
    getopt 'permute' behavior) rather than stopping at the first
    non-option argument, so anything appended after the script name (e.g.
    '--workdir /path') gets fed back into qsub's own option parser and
    rejected as an unrecognized option. We therefore pass all information
    the submitted script needs (working directory, pickled function/kwarg
    paths, etc.) as real job environment variables via qsub's '-v' option
    instead, and route execution through the generic
    `runscripts/pbs_entry_point` wrapper, which reads them back out of the
    environment and translates them into the CLI arguments the underlying
    `submit`/`run` scripts expect.

TODO
    Consider adding support for TORQUE-style 'qstat -f' parsing if/when a
    TORQUE-based child class is required.
"""
import os
import re
import sys
import time
import subprocess
import numpy as np

from seisflows import logger, ROOT_DIR
from seisflows.system.cluster import Cluster
from seisflows.tools import msg
from seisflows.tools.config import pickle_function_list, copy_file


class Pbs(Cluster):
    """
    System PBS
    ----------
    Interface for submitting and monitoring jobs on HPC systems running the
    Portable Batch System (PBS), e.g., PBS Professional, OpenPBS or TORQUE.

    Parameters
    ----------
    :type queue: str
    :param queue: Name of the PBS queue (or routing queue) that jobs are
        submitted to, set with the qsub '-q' option. Must be overwritten by
        child classes with a machine-specific set of available queues (see
        `_queues`).
    :type group: str
    :param group: project/account code used to charge core-hours, set with
        the qsub '-W group_list=<group>' option (on some PBS Professional
        systems this is referred to as the 'project code'). Required by
        most allocation-tracking PBS systems, but left optional here since
        not all PBS deployments require it.
    :type pbs_args: str
    :param pbs_args: Any (optional) additional PBS arguments that will be
        passed to the qsub calls. Should be in the form:
        '-l key1=value1 -l key2=value2'

    Paths
    -----
    ***
    """
    __doc__ = Cluster.__doc__ + __doc__

    # PBS-safe generic entry point wrapper, see module docstring note above.
    # Overwrites `Cluster.submit_workflow`/`Cluster.run_functions`, which
    # point directly at 'submit'/'run' -- unsafe for PBS since neither
    # `Pbs.submit()` nor `Pbs.run()` append CLI arguments after this
    # executable, relying instead on job environment variables (`-v`)
    submit_workflow = os.path.join(ROOT_DIR, "system", "runscripts",
                                   "pbs_entry_point")
    run_functions = os.path.join(ROOT_DIR, "system", "runscripts",
                                 "pbs_entry_point")

    def __init__(self, ntask_max=100, queue=None, group=None, pbs_args="",
                 **kwargs):
        """
        PBS-specific setup parameters

        :type ntask_max: int
        :param ntask_max: set the maximum number of simultaneously running
            array job (sub)processes that are submitted to a cluster at one
            time.
        """
        super().__init__(**kwargs)

        # Overwrite the existing 'mpiexec'
        if self.mpiexec is None:
            self.mpiexec = "mpiexec"
        self.ntask_max = ntask_max
        self.queue = queue
        self.group = group
        self.pbs_args = pbs_args

        # Must be overwritten by child class
        self.submit_to = None
        self._queues = {}

        # Define PBS-dependent job states used for monitoring the queue.
        # 'COMPLETED'/'FAILED' are not native PBS states -- they are
        # synthesized by `query_job_states` from a finished job's
        # 'Exit_status', since PBS's own terminal state (e.g., 'F') does not
        # by itself distinguish success from failure.
        self._completed_states = ["COMPLETED"]
        self._failed_states = ["FAILED"]
        self._pending_states = ["Q", "R", "H", "QUEUED", "RUNNING", "HELD",
                                "BEGUN", "EXITING", "TRANSIT", "WAITING",
                                "SUSPEND"]

    def check(self):
        """
        Checks parameters and paths
        """
        super().check()

        assert(self._queues), (
            f"PBS system child classes require defining `_queues`, a "
            f"lookup of queue name to node_size (cores-per-node) inherent "
            f"to the compute system.")

        assert(self.queue in self._queues), \
            f"PBS `queue` must match one of {list(self._queues.keys())}"

        assert(self.submit_to in self._queues), \
            f"PBS `submit_to` queue must match one of " \
            f"{list(self._queues.keys())}"

    @property
    def nodes(self):
        """Defines the number of nodes which is derived from system node
        size"""
        _nodes = np.ceil(self.nproc / float(self.node_size))
        return int(_nodes)

    @property
    def node_size(self):
        """Defines the node size (cores-per-node) of a given queue. This is a
        hard set number defined by the system architecture"""
        return self._queues[self.queue]

    @property
    def submit_call_header(self):
        """
        The submit call defines the PBS header which is used to submit a
        workflow task list to the system. It is usually dictated by the
        system's required parameters, such as project codes and queues.
        Submit calls are modified and called by the `submit` function.

        :rtype: str
        :return: the system-dependent portion of a submit call
        """
        walltime = self._fmt_walltime(self.walltime)

        _call = " ".join([
            f"qsub",
            f"{self.pbs_args or ''}",
            f"-N {self.title}",
            f"-o {self.path.output_log}",
            f"-j oe",
            f"-q {self.submit_to}",
            f"-l select=1:mpiprocs=1",
            f"-l walltime={walltime}",
        ])
        if self.group:
            _call += f" -W group_list={self.group}"

        return _call

    def run_call(self, executable="", single=False, array=None,
                 tasktime=None, variables=""):
        """
        The run call defines the PBS directives which are used to run tasks
        during an executing workflow. Like the submit call its arguments are
        dictated by the given system. Run calls are modified and called by
        the `run` function.

        .. note::
            Unlike some SLURM-like schedulers, single-task jobs are
            submitted WITHOUT the '-J' array directive (some PBS
            implementations reject single-element array ranges like
            '0-0'). In that case `SEISFLOWS_TASKID` is explicitly set to 0
            via `variables` (`-v`). For actual array jobs, the running task
            instead recovers its ID from the PBS-assigned 'PBS_ARRAY_INDEX'
            environment variable (see `tools.config.ENV_VARIABLES`).

        :type executable: str
        :param exectuable: the actual exectuable to run within the PBS
            directive. Something like './script.py'. No CLI arguments
            should be appended here -- see module docstring note; use
            `variables` instead
        :type array: str
        :param array: overwrite the `array` variable to run specific jobs. If
            not provided, then we will run jobs 0-{ntask}%{ntask_max}. Jobs
            should be submitted in the format of a PBS array string, e.g.,
            0-79%15 or 2-4
        :type single: bool
        :param single: flag to get a run call that is meant to be run on the
            mainsolver (ntask==1), or run for all jobs (ntask times). Examples
            of single process runs include smoothing, and kernel combination
        :type tasktime: float
        :param tasktime: Custom tasktime in units minutes for running the
            given functions. If not given, defaults to the System variable
            `tasktime`
        :type variables: str
        :param variables: comma-separated 'VAR=val' string of job
            environment variables to export via qsub's '-v' option
        :rtype: str
        :return: the system-dependent portion of a run call
        """
        array = array or self.task_ids(single=single)
        tasktime = tasktime or self.tasktime
        walltime = self._fmt_walltime(tasktime)

        use_array = not (single or self.ntask == 1)
        mpiprocs = min(self.nproc, self.node_size)

        _call_list = [
            f"qsub",
            f"{self.pbs_args or ''}",
            f"-N {self.title}",
            f"-j oe",
            f"-o {self.path.log_files}{os.sep}",
            f"-q {self.queue}",
            f"-l select={self.nodes}:mpiprocs={mpiprocs}",
            f"-l walltime={walltime}",
        ]
        if self.group:
            _call_list.append(f"-W group_list={self.group}")
        if use_array:
            # PBS array jobs must be submitted rerunnable ('-r y'), matching
            # e.g. the Miyabi User's Guide's own array job example ('qsub -r
            # y -J ...'). Without this, sites that default the 'Rerunable'
            # job attribute to False will reject array submissions with
            # "cannot submit non-rerunable Array Job"
            _call_list.append("-r y")
            _call_list.append(f"-J {array}")
        if variables:
            _call_list.append(f"-v {variables}")
        _call_list.append(f"{executable}")  # <-- script to run, no CLI args

        return " ".join(_call_list)

    @staticmethod
    def _fmt_walltime(minutes):
        """
        PBS '-l walltime=' and '-l elapse=' resource requests expect a
        string in the format '[[hour:]minute:]second'. We convert from
        SeisFlows' native units of minutes to a 'HH:MM:SS' string.

        :type minutes: float
        :param minutes: walltime/tasktime in units of minutes
        :rtype: str
        :return: formatted time string 'HH:MM:SS'
        """
        return "{:02d}:{:02d}:00".format(*divmod(round(minutes), 60))

    @staticmethod
    def _stdout_to_job_id(stdout):
        """
        The stdout message after a `qsub` job is submitted is simply the
        assigned Job ID, e.g., '123456.opbs', or for array jobs,
        '123456[].opbs'.

        :type stdout: str
        :param stdout: standard qsub response after submitting a job
        :rtype: str
        :return: the matching Job ID
        :raises SystemExit: if the job id does not match the expected format
        """
        job_id = str(stdout).strip().split("\n")[-1].strip()
        if not re.match(r"^\d+(\[])?\.\S+$", job_id):
            logger.critical(f"parsed job id '{job_id}' does not match the "
                            f"expected PBS job id format "
                            f"'<sequence>[.<index>].<server>', please check "
                            f"that function `system._stdout_to_job_id()` is "
                            f"set correctly")
            sys.exit(-1)

        return job_id

    def _extra_qsub_variables(self):
        """
        Hook for child classes to inject additional job environment
        variables (beyond the ones `Pbs.submit()`/`Pbs.run()` already set)
        into the '-v' option of a qsub call, e.g. Miyabi uses this to pass
        `SEISFLOWS_CONDA_ENV` to `runscripts/conda_activate-miyabi`.

        :rtype: str
        :return: comma-separated 'VAR=val' string, or empty string
        """
        return ""

    @property
    def _runscripts_dir(self):
        """
        Absolute path to `system/runscripts/`, passed to submitted jobs as
        `SEISFLOWS_RUNSCRIPTS_DIR` (see `submit`/`run`).

        .. note::
            qsub COPIES a submitted script into its own spool directory
            (e.g. '/var/spool/pbs/mom_priv/jobs/') before executing it, so
            a running script cannot reliably locate its sibling scripts by
            self-referencing its own path (e.g. bash's
            `dirname "${BASH_SOURCE[0]}"`) -- that would resolve to the
            spool directory, not this one. We therefore pass this path
            explicitly instead.

        :rtype: str
        :return: absolute path to the runscripts directory
        """
        return os.path.join(ROOT_DIR, "system", "runscripts")

    def submit(self, workdir=None, parameter_file="parameters.yaml"):
        """
        Submits the main workflow job as a separate job submitted directly
        to the system that is running the master job.

        .. note::
            Overwrites `Cluster.submit()`. See module docstring note --
            `--workdir`/`--parameter_file` are passed via qsub's '-v'
            option (as `SEISFLOWS_WORKDIR`/`SEISFLOWS_PARAMETER_FILE`)
            rather than as trailing CLI arguments to `submit_workflow`.

        :type workdir: str
        :param workdir: path to the current working directory
        :type parameter_file: str
        :param parameter_file: parameter file name used to instantiate the
            SeisFlows package
        """
        # Copy log files if present to avoid overwriting
        for src in [self.path.output_log, self.path.par_file]:
            if os.path.exists(src) and os.path.exists(self.path.log_files):
                copy_file(src, copy_to=self.path.log_files)

        workdir = workdir or self.path.workdir
        variables = (f"SEISFLOWS_ENTRY_POINT=submit,"
                    f"SEISFLOWS_RUNSCRIPTS_DIR={self._runscripts_dir},"
                    f"SEISFLOWS_WORKDIR={workdir},"
                    f"SEISFLOWS_PARAMETER_FILE={parameter_file}")
        extra = self._extra_qsub_variables()
        if extra:
            variables += f",{extra}"

        submit_call = " ".join([
            self.submit_call_header,
            f"-v {variables}",
            f"{self.submit_workflow}",
        ])
        logger.debug(submit_call)
        try:
            subprocess.run(submit_call, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.critical(f"SeisFlows master job has failed with: {e}")
            sys.exit(-1)

    def run(self, funcs, single=False, tasktime=None, array=None,
            _attempts=0, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fashion on a PBS
        cluster. Executes the list of functions (`funcs`) NTASK times with
        each task occupying NPROC cores.

        .. note::

            Completely overwrites the `Cluster.run()` command

        :type funcs: list of methods
        :param funcs: a list of functions that should be run in order. All
            kwargs passed to run() will be passed into the functions.
        :type single: bool
        :param single: run a single-process, non-parallel task, such as
            smoothing the gradient, which only needs to be run by once.
            This will change how the job array and the number of tasks is
            defined, such that the job is submitted as a single-core job to
            the system.
        :type tasktime: float
        :param tasktime: Custom tasktime in units minutes for running the
            given functions `funcs`. If not given, defaults to the System
            variable `tasktime`. If tasks exceed the given `tasktime`, the
            program will exit
        :type array: str
        :param array: overwrite the `array` variable to run specific jobs. If
            not provided, then we will run jobs 0-{ntask}%{ntask_max}. Jobs
            should be submitted in the format of a PBS array string, e.g.,
            0-79%15
        :type _attempts: int
        :param _attempts: a recursive counter for failed job runs that allows
            the `run` function to re-attempt failed jobs up to `rerun` number
            of times
        """
        logger.info(f"running functions {[_.__name__ for _ in funcs]}")

        # Condense the functions that we want to run on system to a pickle
        # file
        funcs_fid, kwargs_fid = pickle_function_list(
                funcs, path=self.path.scratch, verbose=self.verbose,
                level=self.log_level, **kwargs
                )

        # Pass everything the running task needs via job environment
        # variables ('-v') rather than as CLI arguments -- see module
        # docstring note. `SEISFLOWS_TASKID` is only needed for non-array
        # (single-task) jobs; array (sub)jobs recover it from PBS's own
        # 'PBS_ARRAY_INDEX' environment variable instead. `SEISFLOWS_WORKDIR`
        # is used by `runscripts/pbs_entry_point` to `cd` into the working
        # directory before running -- unlike SLURM's sbatch, PBS jobs do NOT
        # default to the submission directory (PBS User's Guides instead
        # recommend an explicit 'cd ${PBS_O_WORKDIR}' in job scripts)
        variables = (f"SEISFLOWS_ENTRY_POINT=run,"
                    f"SEISFLOWS_RUNSCRIPTS_DIR={self._runscripts_dir},"
                    f"SEISFLOWS_WORKDIR={self.path.workdir},"
                    f"SEISFLOWS_FUNCS={funcs_fid},"
                    f"SEISFLOWS_KWARGS={kwargs_fid}")
        if single or self.ntask == 1:
            variables += ",SEISFLOWS_TASKID=0"
        if self.environs:
            variables += f",{self.environs}"
        extra = self._extra_qsub_variables()
        if extra:
            variables += f",{extra}"

        # Get the run call that will be submitted to the system via
        # subprocess
        run_call = self.run_call(executable=f"{self.run_functions}",
                                 tasktime=tasktime, array=array,
                                 single=single, variables=variables)
        logger.debug(run_call)

        # RUN the job by submitting the qsub directive to system
        try:
            stdout = subprocess.run(run_call, stdout=subprocess.PIPE,
                                    text=True, shell=True,
                                    check=True).stdout
        except subprocess.CalledProcessError as e:
            logger.critical(msg.cli(
                f"SeisFlows task submission has failed with: {e}",
                items=[f"RUN CALL: {run_call}"],
                header="pbs run error", border="=")
                )
            sys.exit(-1)

        job_id = self._stdout_to_job_id(stdout)

        # Monitor the job queue until all jobs have finished
        status = self.monitor_job_status(job_id)

        # Failed job can either try to re-run, or simply exit the program
        if status == -1:
            # Failure recovery mechanism. Determine which jobs did not
            # complete and attempt to rerun them up to a certain amount of
            # times
            if self.rerun and _attempts < self.rerun:
                jobs, states = self.query_job_states(job_id, sort=True)
                # Assumes that the (sorted) job index matches the task IDs
                failed_array = []
                for i, (job, state) in enumerate(zip(jobs, states)):
                    if state in self._failed_states:
                        failed_array.append(i)
                array_str = ",".join([str(_) for _ in failed_array])

                logger.info(f"attempt {_attempts+1}/{self.rerun} rerun "
                            f"{len(failed_array)} failed jobs")
                logger.debug(f"task ids to rerun: {array_str}")

                # Recursively 'run' the functions but only with the failed
                # jobs, assuming that all other jobs completed nominally
                self.run(funcs=funcs, single=single, tasktime=tasktime,
                         array=array_str, _attempts=_attempts + 1, **kwargs)
            else:
                logger.critical(
                    msg.cli(f"Stopping workflow. Please check logs for "
                            f"details",
                            items=[f"TASKS:   {[_.__name__ for _ in funcs]}",
                                   f"JOB_ID: {job_id}",
                                   f"QSUB:  {run_call}"],
                            header="pbs run error", border="=")
                )
                sys.exit(-1)
        # Completed job will end the 'run' function
        else:
            logger.info(f"task {job_id} finished successfully")
            # Wait for all processes to finish and write to disk (if they
            # do). Moving on too quickly may result in required files not
            # being available
            time.sleep(5)

    def task_ids(self, single=False):
        """
        Overwrite `system.workstation.task_ids` to get PBS specific array
        configurations which are passed as strings for the '-J' PBS
        directive, rather than lists which is how `system.workstation`
        handles this.

        .. note::
            Unlike SLURM, PBS's '-J' array directive only accepts a single
            contiguous range with an optional step
            (`<start>-<end>[:step]`), not an arbitrary comma-separated list
            of ranges. If `system.array` is set by the User, it is assumed
            to already be given in PBS-compatible format.

        :type single: bool
        :param single: If we only want to run a single process, this will
            default to TaskID == 0
        :rtype: str
        :return: string formatter of Task IDs to be used by the `run`
            function via the `run_call`
        """
        if single:
            task_ids = "0-0"
        else:
            if self.array is not None:
                task_ids = self.array
            else:
                # e.g., 0-79%15 means 80 jobs, 15 at a time
                task_ids = f"0-{self.ntask - 1}%{self.ntask_max}"

        return task_ids

    def query_job_states(self, job_id, sort=False):
        """
        Overwrites `system.cluster.Cluster.query_job_states`

        Queries completion status of a job (or job array) using `qstat`.
        Because PBS's `qstat` only lists jobs that are queued/running/held
        (not finished ones), and because a finished job's terminal STATUS
        does not by itself distinguish success from failure, this function:

        1) Queries currently active (queued/running/held) jobs with `qstat`
        2) If none are found (i.e., the job has left the active queue),
           queries the finished-job history with `qstat -H`
        3) For any jobs found to be finished, queries `tracejob` to
           determine each job's 'Exit_status' (0 == success), and
           synthesizes 'COMPLETED'/'FAILED' states from this

        :type job_id: str
        :param job_id: main job id (or array job id) to query, returned
            from the subprocess.run that ran the job(s)
        :type sort: bool
        :param sort: sort by (sub)job array index. Defaults to False because
            currently running jobs may return job numbers that cannot be
            sorted, e.g., '1[].opbs' while queued. We only use sort when
            recovering from job failure because then we are assured that
            all (sub)jobs have run.
        :rtype: (list, list)
        :return: (job ids, corresponding job states). Returns (None, None)
            if no information can be found for `job_id` (e.g., jobs have not
            yet initialized on system)
        """
        job_ids, job_states = self._query_active_job_states(job_id)

        # Jobs are still queued/running/held, nothing more to do
        if job_ids:
            if sort:
                job_ids, job_states = self._sort_by_array_index(job_ids,
                                                                 job_states)
            return job_ids, job_states

        # No active jobs found, check whether the job(s) have finished
        finished_ids = self._query_finished_job_ids(job_id)
        if not finished_ids:
            return None, None

        # Determine success/failure of each finished (sub)job via
        # `tracejob`'s reported 'Exit_status' (0 == success)
        for jid in finished_ids:
            cmd = f"tracejob -n 1 {jid}"
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    shell=True)
            match = re.search(r"Exit_status=(-?\d+)", result.stdout)
            job_ids.append(jid)
            if match and int(match.group(1)) == 0:
                job_states.append("COMPLETED")
            else:
                job_states.append("FAILED")

        if sort:
            job_ids, job_states = self._sort_by_array_index(job_ids,
                                                             job_states)

        return job_ids, job_states

    @staticmethod
    def _parse_qstat_stdout(stdout):
        """
        Parse the JOB_ID and STATUS columns out of `qstat`-formatted stdout
        (columns: JOB_ID, JOB_NAME, STATUS, ...)

        :type stdout: str
        :param stdout: stdout from a `qstat` call
        :rtype: (list, list)
        :return: (list of job ids, list of corresponding job states)
        """
        job_ids, job_states = [], []
        for line in str(stdout).strip().splitlines():
            line = line.strip()
            # Skip header/banner lines, e.g. 'JOB_ID JOB_NAME STATUS ...'
            # or 'Miyabi scheduled stop time: ...'
            if not line or not line[0].isdigit():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            job_ids.append(parts[0])
            job_states.append(parts[2].upper())

        return job_ids, job_states

    def _query_active_job_states(self, job_id):
        """
        Query currently queued/running/held (sub)jobs with `qstat -t`

        :type job_id: str
        :param job_id: main job id (or array job id) to query
        :rtype: (list, list)
        :return: (list of job ids, list of corresponding job states)
        """
        cmd = f"qstat -t {job_id}"
        result = subprocess.run(cmd, capture_output=True, text=True,
                                shell=True)
        return self._parse_qstat_stdout(result.stdout)

    def _query_finished_job_ids(self, job_id):
        """
        Query finished (sub)job ids with `qstat -H`

        :type job_id: str
        :param job_id: main job id (or array job id) to query
        :rtype: list
        :return: list of finished (sub)job ids
        """
        cmd = f"qstat -H -t {job_id}"
        result = subprocess.run(cmd, capture_output=True, text=True,
                                shell=True)
        finished_ids, _ = self._parse_qstat_stdout(result.stdout)

        return finished_ids

    @staticmethod
    def _sort_by_array_index(job_ids, job_states):
        """
        Sort (sub)job ids/states by their array index, e.g., so that
        '123456[9].opbs' sorts after '123456[2].opbs'. (Sub)job ids without
        an array index (i.e., non-array jobs) sort first.

        :type job_ids: list
        :param job_ids: list of job ids to sort
        :type job_states: list
        :param job_states: list of job states corresponding to `job_ids`
        :rtype: (tuple, tuple)
        :return: (sorted job ids, correspondingly sorted job states)
        """
        def _idx(jid):
            match = re.search(r"\[(\d+)]", jid)
            return int(match.group(1)) if match else -1

        job_ids, job_states = zip(
            *sorted(zip(job_ids, job_states), key=lambda x: _idx(x[0]))
        )

        return job_ids, job_states
