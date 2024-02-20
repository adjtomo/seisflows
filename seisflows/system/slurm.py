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
from seisflows import logger
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
        self._walltime = str(timedelta(minutes=self.walltime))
    
        # Define SLURM-dependent job states used for monitoring queue
        # Available job states are listed here: 
        #                   https://slurm.schedmd.com/sacct.html
        self._completed_states = ["COMPLETED"]
        self._failed_states = ["TIMEOUT", "FAILED", "NODE_FAIL", 
                               "OUT_OF_MEMORY", "CANCELLED"]
        self._pending_states = ["PENDING", "RUNNING"]

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

        assert("--parsable" in self.run_call()), (
            f"System `run_call` requires SBATCH argument '--parsable' "
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

    def run_call(self, executable="", single=False, array=None, tasktime=None):
        """
        The run call defines the SBATCH call which is used to run tasks during
        an executing workflow. Like the submit call its arguments are dictated
        by the given system. Run calls are modified and called by the `run`
        function

        :type executable: str
        :param exectuable: the actual exectuable to run within the SBATCH 
            directive. Something like './script.py'
        :type array: str
        :param array: overwrite the `array` variable to run specific jobs. If
            not provided, then we will run jobs 0-{ntask}%{ntask_max}. Jobs 
            should be submitted in the format of a SLURM array string, 
            something like: 0,1,3,5 or 2-4,8-22
        :type single: bool
        :param single: flag to get a run call that is meant to be run on the
            mainsolver (ntask==1), or run for all jobs (ntask times). Examples
            of single process runs include smoothing, and kernel combination
        :rtype: str
        :return: the system-dependent portion of a run call
        """
        array = array or self.task_ids(single=single)  # get job array str
        if tasktime is None:
            # 0 is there just for initialization. `tasktime` will superceded
            tasktime = str(timedelta(minutes=0 or self.tasktime))
        else:
            tasktime = str(timedelta(minutes=tasktime))

        # Determine if this is a single-process or array job
        if single:
            ntasks = 1
            env = "SEISFLOWS_TASKID=0," 
        else:
            ntasks = self.nproc
            env = ""
        
        _call = " ".join([
             f"sbatch",
             f"{self.slurm_args or ''}",
             f"--job-name={self.title}",
             f"--nodes={self.nodes}",
             f"--ntasks-per-node={self.node_size:d}",
             f"--ntasks={ntasks:d}",
             f"--time={tasktime}",
             f"--output={os.path.join(self.path.log_files, '%A_%a')}",
             f"--array={array}",
             f"--parsable",
             f"{executable}",  # <-- The actual script/program to run goes here
             f"--environment {env}{self.environs or ''}"
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
        except ValueError as e:
            logger.critical(f"parsed job id '{job_id}' does not evaluate as an "
                            f"integer, please check that function "
                            f"`system._stdout_to_job_id()` is set correctly")
            sys.exit(-1)

        return job_id

    def run(self, funcs, single=False, tasktime=None, array=None, _attempts=0, 
            **kwargs):
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
        :type array: str
        :param array: overwrite the `array` variable to run specific jobs. If
            not provided, then we will run jobs 0-{ntask}%{ntask_max}. Jobs 
            should be submitted in the format of a SLURM array string, 
            something like: 0,1,3,5 or 2-4,8-22
        :type _attempts: int
        :param _attempts: a recursive counter for failed job runs that allows 
            the `run` function to re-attempt failed jobs up to `rerun` number
            of times
        """
        logger.info(f"running functions {[_.__name__ for _ in funcs]}")

        # Condense the functions that we want to run on system to a pickle file
        funcs_fid, kwargs_fid = pickle_function_list(
                funcs, path=self.path.scratch, verbose=self.verbose,
                level=self.log_level, **kwargs
                )

        # Get the run call that will be submitted to the system via subprocess
        run_call = self.run_call(executable=f"{self.run_functions} "
                                            f"--funcs {funcs_fid} "
                                            f"--kwargs {kwargs_fid}",
                                 tasktime=tasktime, array=array, single=single
                                 )
        logger.debug(run_call)

        # RUN the job by submitting the sbatch directive to system
        stdout = subprocess.run(run_call, stdout=subprocess.PIPE,
                                text=True, shell=True).stdout
        job_id = self._stdout_to_job_id(stdout)

        # Monitor the job queue until all jobs have finished
        status = self.monitor_job_status(job_id)

        # Failed job can either try to re-run, or simply exit the program
        if status == -1:  
            # Failure recovery mechanism. Determine which jobs did not complete
            # and attempt to rerun them up to a certain amount of times
            if self.rerun and _attempts < self.rerun:
                jobs, states = self.query_job_states(job_id, sort=True)
                # Here we are assuming that the job index matches the task IDs
                failed_array = []
                for i, (job, state) in enumerate(zip(jobs, states)):
                    if state in self._failed_states:
                        failed_array.append(i)
                array_str = ",".join([str(_) for _ in failed_array])

                logger.info(f"attempt {_attempts+1}/{self.rerun} rerun "
                            f"{len(failed_array)} failed jobs")
                logger.debug(f"task ids to rerun: {array_str}")

                # Recursively 'run' the functions but only with the failed jobs,
                # assuming that all other jobs completed nominally
                self.run(funcs=funcs, single=single, tasktime=tasktime,
                         array=array_str, _attempts=_attempts+1, **kwargs)
            else:
                logger.critical(
                    msg.cli(f"Stopping workflow. Please check logs for details",
                            items=[f"TASKS:   {[_.__name__ for _ in funcs]}",
                                   f"JOB_ID: {job_id}",
                                   f"SBATCH:  {run_call}"],
                            header="slurm run error", border="=")
                )
                sys.exit(-1)
        # Completed job will end the 'run' function
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
            via the `run_call`
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

    def query_job_states(self, job_id, sort=False):
        """
        Overwrites `system.cluster.Cluster.query_job_states`

        Queries completion status of an array job by running the SLURM `sacct`

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
        :param job_id: main job id to query, returned from the subprocess.run 
            that ran the jobs
        :type sort: bool
        :param sort: sort by job ids or job array ids. Defaults to False because
            currently running jobs may return job numbers that cannot be sorted
            e.g., 1_0, 1_1, 1_[2-5]. We only use sort when recovering from job 
            failure because then we are assured that all jobs have run.
        :rtype: (list, list)
        :return: (job ids, corresponding job states). Returns (None, None) if
            `sacct` does not return a useful stdout (e.g., jobs have not
            yet initialized on system)
        """
        job_ids, job_states = [], []
        cmd = f"sacct -nLX -o jobid,state -j {job_id}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        stdout = result.stdout

        # If no return, return None so calling function knows something is wrong
        if not stdout:
            return None, None

        # Return the job numbers and respective states for the given job ID
        for job_line in str(stdout).strip().split("\n"):
            if not job_line:
                continue
            job_id, job_state = job_line.split()
            job_ids.append(job_id)
            job_states.append(job_state)
    
        if sort:
            # Sort by job ids because we assume that logically job numbers are 
            # in a numerically ascending order i.e., 1,2,3 
            job_ids, job_states = zip(*sorted(zip(job_ids, job_states)))

            # Sort array jobs because normal 'sorted' function doesn't work when
            # strings are hyphenated (e.g., 1_0, 1_1, 1_2)
            # https://stackoverflow.com/questions/20862968/\
            #            numbers-with-hyphens-or-strings-of-numbers-with-hyphens
            job_ids, job_states = zip(
                *sorted(zip(job_ids, job_states), 
                        key=lambda x: [int(y) for y in x[0].split('_')])
                        )

        return job_ids, job_states

