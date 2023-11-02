#!/usr/bin/env python3
"""
The Simple Linux Utility for Resource Management (SLURM) is a commonly used
workload manager on many high performance computers / clusters. The Slurm
system class provides generalized utilites for interacting with Slurm systems.

Useful commands for figuring out system-specific required parameters
    $ sinfo --Node --long  # Determine the cores-per-node for partitions

.. note::
    The main development system for SeisFlows used SLURM. Therefore the other
    system supers will not be up to date until access to those systems are
    granted. This rosetta stone, for converting from SLURM to other workload
    management tools will be useful: https://slurm.schedmd.com/rosetta.pdf

.. note::
   SLURM systems expect walltime/tasktime in format: "minutes", 
   "minutes:seconds", "hours:minutes:seconds". SeisFlows uses the latter
   and converts task and walltimes from input of minutes to a time string.

TODO
    Create 'slurm_singulairty', a child class for singularity-based runs which
    loads and runs programs through singularity, OR add a parameter options
    which will change the run and/or submit calls
"""
import os
import sys
import numpy as np
import time
import subprocess

from datetime import timedelta
from seisflows import ROOT_DIR, logger
from seisflows.system.cluster import Cluster
from seisflows.tools import msg
from seisflows.tools.config import pickle_function_list


class Slurm(Cluster):
    """
    System Slurm
    ------------
    Interface for submitting and monitoring jobs on HPC systems running the 
    Simple Linux Utility for Resource Management (SLURM) workload manager.

    Parameters
    ----------
    :type slurm_args: str
    :param slurm_args: Any (optional) additional SLURM arguments that will
        be passed to the SBATCH scripts. Should be in the form:
        '--key1=value1 --key2=value2"

    Paths
    -----
    ***
    """
    __doc__ = Cluster.__doc__ + __doc__

    def __init__(self, ntask_max=100, slurm_args="", **kwargs):
        """
        Slurm-specific setup parameters

        :type ntask_max: int
        :param ntask_max: set the maximum number of simultaneously running array
            job processes that are submitted to a cluster at one time.
        """
        super().__init__(**kwargs)

        # Overwrite the existing 'mpiexec'
        if self.mpiexec is None:
            self.mpiexec = "srun -u"
        self.ntask_max = ntask_max
        self.slurm_args = slurm_args

        # Must be overwritten by child class
        self.partition = None
        self.submit_to = None
        self._partitions = {}

        # Convert walltime and tasktime to datetime str 'H:MM:SS'
        self._tasktime = None
        self._walltime = str(timedelta(minutes=self.walltime))
    
        # Define SLURM-dependent job states used for monitoring queue
        # Available job states are listed here: 
        #                   https://slurm.schedmd.com/sacct.html
        self._completed_states = ["COMPLETED"]
        self._failed_states = ["TIMEOUT", "FAILED", "NODE_FAIL", 
                               "OUT_OF_MEMORY", "CANCELLED"]
        self._pending_states = ["PENDING"]

    def check(self):
        """
        Checks parameters and paths
        """
        super().check()

        assert(self.node_size is not None), (
            f"Slurm system child classes require defining the node_size or "
            f"the number of cores per node inherent to the compute system.")

        assert(self.partition in self._partitions), \
            f"Cluster partition name must match {self._partitions}"

        assert(self.submit_to in self._partitions), \
            f"Submission partition name must match {self._partitions}"

        assert("--parsable" in self.run_call_header), (
            f"System `run_call_header` requires SBATCH argument '--parsable' "
            f"which is required to keep STDOUT formatted correctly when "
            f"submitting jobs to the system."
        )

    @property
    def nodes(self):
        """Defines the number of nodes which is derived from system node size"""
        _nodes = np.ceil(self.nproc / float(self.node_size))
        _nodes = _nodes.astype(int)
        return _nodes

    @property
    def node_size(self):
        """Defines the node size of a given cluster partition. This is a hard
        set number defined by the system architecture"""
        return self._partitions[self.partition]

    @property
    def submit_call_header(self):
        """
        The submit call defines the SBATCH header which is used to submit a
        workflow task list to the system. It is usually dictated by the
        system's required parameters, such as account names and partitions.
        Submit calls are modified and called by the `submit` function.

        :rtype: str
        :return: the system-dependent portion of a submit call
        """
        _call = " ".join([
            f"sbatch",
            f"{self.slurm_args or ''}",
            f"--job-name={self.title}",
            f"--output={self.path.output_log}",
            f"--error={self.path.output_log}",
            f"--ntasks-per-node={self.node_size}",
            f"--nodes=1",
            f"--time={self._walltime}"
        ])
        return _call

    @property
    def run_call_header(self):
        """
        The run call defines the SBATCH header which is used to run tasks during
        an executing workflow. Like the submit call its arguments are dictated
        by the given system. Run calls are modified and called by the `run`
        function

        :rtype: str
        :return: the system-dependent portion of a run call
        """
        _call = " ".join([
             f"sbatch",
             f"{self.slurm_args or ''}",
             f"--job-name={self.title}",
             f"--nodes={self.nodes}",
             f"--ntasks-per-node={self.node_size:d}",
             f"--ntasks={self.nproc:d}",
             f"--time={self._tasktime}",
             f"--output={os.path.join(self.path.log_files, '%A_%a')}",
             f"--array={self.task_ids()}",
             f"--parsable"
        ])
        return _call

    @staticmethod
    def _stdout_to_job_id(stdout):
        """
        The stdout message after an SBATCH job is submitted, from which we get
        the job number, differs between systems, allow this to vary

        .. note:: Examples
            1) standard example: Submitted batch job 4738244
            2) (1) with '--parsable' flag: 4738244
            3) federated cluster: Submitted batch job 4738244; Maui
            4) (3) with '--parsable' flag: 4738244; Maui
        
        This function deals with cases (2) and (4). Other systems that have more 
        complicated stdout messages will need to overwrite this function

        :type stdout: str
        :param stdout: standard SBATCH response after submitting a job with the
            '--parsable' flag
        :rtype: str
        :return: a matching job ID. We convert str->int->str to ensure that
            the job id is an integer value (which it must be)
        :raises SystemExit: if the job id does not evaluate as an integer
        """
        job_id = str(stdout).split(";")[0].strip()
        try:
            int(job_id)
        except ValueError:
            logger.critical(f"parsed job id '{job_id}' does not evaluate as an "
                            f"integer, please check that function "
                            f"`system._stdout_to_job_id()` is set correctly")
            sys.exit(-1)

        return job_id

    def run(self, funcs, single=False, tasktime=None, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on a SLURM
        cluster. Executes the list of functions (`funcs`) NTASK times with each
        task occupying NPROC cores.

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
        :param tasktime: Custom tasktime in units minutes for running the given
            functions `funcs`. If not given, defaults to the System variable
            `tasktime`. If tasks exceed the given `tasktime`, the program will 
            exit
        """
        # Allow custom tasktime, else default to System variable
        if tasktime is None:
            tasktime = self.tasktime
        self._tasktime = str(timedelta(minutes=tasktime))

        funcs_fid, kwargs_fid = pickle_function_list(
                funcs, path=self.path.scratch, verbose=self.verbose,
                level=self.log_level, **kwargs
                )

        if single:
            logger.info(f"running functions {[_.__name__ for _ in funcs]} on "
                        f"system 1 time")
        else:
            logger.info(f"running functions {[_.__name__ for _ in funcs]} on "
                        f"system {self.ntask} times")

        # Default sbatch command line input, can be overloaded by subclasses
        # Copy-paste this default run_call and adjust accordingly for subclass
        run_call = " ".join([
            f"{self.run_call_header}",
            f"{self.run_functions}",
            f"--funcs {funcs_fid}",
            f"--kwargs {kwargs_fid}",
            f"--environment {self.environs or ''}"
        ])

        # Single-process jobs simply need to replace a few sbatch arguments.
        # Do it AFTER `run_call` has been defined so that subclasses submitting
        # custom run calls can still benefit from this
        if single:
            logger.debug("replacing parts of sbatch run call for single "
                         "process job")
            run_call = self._modify_run_call_single_proc(run_call)

        logger.debug(run_call)

        # RUN the job by submitting the sbatch directive to system
        # get job id (used to monitor job status) from the stdout message
        stdout = subprocess.run(run_call, stdout=subprocess.PIPE,
                                text=True, shell=True).stdout
        job_id = self._stdout_to_job_id(stdout)

        # Monitor the job queue until all jobs have completed, or any one fails
        try:
            status = self.monitor_job_status(job_id)
        except FileNotFoundError:
            logger.critical(f"cannot access job information through 'sacct', "
                            f"waited {TIMEOUT_S}s with no return, please "
                            f"check job scheduler and log messages, or "
                            f"increase timeout constant in `system.slurm`")
            sys.exit(-1)

        # Job has completed
        if status == -1:  # Failed job
            logger.critical(
                msg.cli(f"Stopping workflow. Please check logs for details.",
                        items=[f"TASKS:   {[_.__name__ for _ in funcs]}",
                               f"SBATCH:  {run_call}"],
                        header="slurm run error", border="=")
            )
            sys.exit(-1)
        else:
            logger.info(f"task {job_id} finished successfully")
            # Wait for all processes to finish and write to disk (if they do)
            # Moving on too quickly may result in required files not being avail
            time.sleep(5)

    def task_ids(self, single=False):
        """
        Overwrite `system.workstation.task_ids` to get SLURM specific array
        configurations which are passed as strings for the --array={task_ids()}
        SLURM directive, rather than lists which is how `system.workstation`
        handles this

        Relevant format definition: https://slurm.schedmd.com/job_array.html

        :type single: bool
        :param single: If we only want to run a single process, this is will
            default to TaskID == 0
        :rtype: str
        :return: string formatter of Task IDs to be used by the `run` function
            via the `run_call_header`
        """
        if single:
            task_ids = "0"
        else:
            if self.array is not None:
                task_ids = self.array
            else:
                # e.g., 0-79%15 means 80 jobs, 15 at a time
                task_ids = f"0-{self.ntask-1}%{self.ntask_max}"

        return task_ids

    def query_job_states(self, job_id, timeout_s=300, wait_time_s=30, 
                         _recheck=0):
        """
        Overwrites `system.cluster.Cluster.query_job_states`

        Queries completion status of an array job by running the SLURM cmd 
        `sacct`, 

        .. note::
            The actual command line call wil look something like this
            $ sacct -nLX -o jobid,state -j 441630
            441630_0    PENDING
            441630_1    COMPLETED

        .. note::
            SACCT flag options are described as follows:
            -L: queries all available clusters, not just the cluster that ran 
                the `sacct` call. Used for federated clusters
            -X: supress the .batch and .extern jobnames that are normally 
                returned but don't represent that actual running job

        :type job_id: str
        :param job_id: main job id to query, returned from the subprocess.run that
            ran the jobs
        :type timeout_s: float
        :param timeout_s: length of time [s] to wait for system to respond 
            nominally before rasing a TimeoutError. Sometimes it takes a while
            for job monitoring to start working.
        :type wait_time_s: float
        :param wait_time_s: initial wait time to allow System to initialize 
            jobs in the queue. I.e., how long do we  expect it to take before 
            the job we submitted shows up in the queue.
        :type rechecks: int
        :param rechecks: Used for recursive calling of the function. It can take 
            time for jobs to be initiated on a system, which may result in the 
            stdout of the 'sacct' command to be empty. In this case we wait and 
            call the function again. Rechecks are used to prevent endless loops 
            by putting a stop criteria
        :raises TimeoutError: if 'sacct' does not return any output for ~1 min.
        """
        job_ids, job_states = [], []
        cmd = f"sacct -nLX -o jobid,state -j {job_id}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        stdout = result.stdout

        # Recursively re-check job state incase the job has not been instantiated 
        # in which cause 'stdout' is an empty string
        if not stdout:
            _recheck += 1
            if _recheck > (timeout_s // wait_time_s):
                raise TimeoutError(f"cannot access job ID {job_id}")
            time.sleep(wait_time_s)
            # Recursive call while ticking up `_recheck` as a timeout counter
            query_job_states(job_id, _recheck=_recheck)

        # Return the job numbers and respective states for the given job ID
        for job_line in str(stdout).strip().split("\n"):
            if not job_line:
                continue
            job_id, job_state = job_line.split()
            job_ids.append(job_id)
            job_states.append(job_state)

        return job_ids, job_states

    def _modify_run_call_single_proc(self, run_call):
        """
        Modifies a SLURM SBATCH command to use only 1 processor as a single run
        by replacing the --array and --ntasks options

        :type run_call: str
        :param run_call: The SBATCH command to modify
        :rtype: str
        :return: a modified SBATCH command that should only run on 1 processor
        """
        for part in run_call.split(" "):
            if "--array" in part:
                run_call = run_call.replace(part, "--array=0")
            elif "--ntasks" in part:
                run_call = run_call.replace(part, "--ntasks=1")

        # Append taskid to environment variable, deal with the case where
        # self.par.ENVIRONS is an empty string
        task_id_str = "SEISFLOWS_TASKID=0"
        if not run_call.strip().endswith("--environment"):
            task_id_str = f",{task_id_str}"  # appending to the list of vars

        run_call += task_id_str

        return run_call

