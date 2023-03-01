#!/usr/bin/env python3
"""
Fujitsu brand clusters run their own Job Management System (Fujitsu Technical 
Computing Suite) which we abbreviate as `PJM`. 

PJM is similar to PBS/TORQUE but has its own unique peculiarities. Systems like 
UToyko's Wisteria run the Fujitsu job scheduler
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


# Define bad states defined by SLURM which signifiy failed jobs
BAD_STATES = ["TIMEOUT", "FAILED", "NODE_FAIL", "OUT_OF_MEMORY", "CANCELLED"]


class Pjm(Cluster):
    """
    System Pjm
    ----------
    Interface for submitting and monitoring jobs on HPC systems running the 
    Fujitsu job management system, abbreviated PJM

    Parameters
    ----------

    Paths
    -----
    ***
    """
    __doc__ = Cluster.__doc__ + __doc__

    def __init__(self, group=None, resource_group=None, ntask_max=100, 
                 pjm_args="",  **kwargs):
        """
        Fujitsu-specific setup parameters

        :type ntask_max: int
        :param ntask_max: set the maximum number of simultaneously running array
            job processes that are submitted to a cluster at one time.
        """
        super().__init__(**kwargs)

        # Overwrite the existing 'mpiexec'
        if self.mpiexec is None:
            self.mpiexec = "mpiexec"
        self.ntask_max = ntask_max
        self.group = group
        self.resource_group = resource_group
        self.pjm_args = pjm_args

        # Must be overwritten by child class
        self.partition = None
        self.submit_to = None
        self._partitions = {}

        # Convert walltime and tasktime to datetime str 'H:MM:SS'
        self._tasktime = str(timedelta(minutes=self.tasktime))
        self._walltime = str(timedelta(minutes=self.walltime))

    def check(self):
        """
        Checks parameters and paths
        """
        super().check()

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
        The submit call defines the PJM header which is used to submit a
        workflow task list to the system. It is usually dictated by the
        system's required parameters, such as account names and partitions.
        Submit calls are modified and called by the `submit` function.

        :rtype: str
        :return: the system-dependent portion of a submit call
        """
        _call = " ".join([
            f"pjsub",
            f"{self.pjm_args or ''}",
            f"-L rscgrp={self.resource_group}",  # resource group
            f"-g {self.self.group}",  # project code
            f"-N {self.title}",  # job name
            f"-o {self.path.output_log}",  # write stdout to file
            f"-e {self.path.output_log}",  # write stderr to file
            f"-L elapse={self._walltime}",  # [[hour:]minute:]second
            f"-L node=1",
            f"--mpi proc=1",
        ])
        return _call

    @property
    def run_call_header(self):
        """
        The run call defines the PJM header which is used to run tasks during
        an executing workflow. Like the submit call its arguments are dictated
        by the given system. Run calls are modified and called by the `run`
        function

        TODO 
            array job?
            partition sizes?

        :rtype: str
        :return: the system-dependent portion of a run call
        """
        _call = " ".join([
             f"pjsub",
             f"{self.pjm_args or ''}",
             f"-L rscgrp={self.resource_group}",  # resource group
             f"-g {self.self.group}",  # project code
             f"-N {self.title}",  # job name
             f"-o {self.path.output_log}",  # write stdout to file
             f"-e {self.path.output_log}",  # write stderr to file
             f"-L elapse={self._walltime}",  # [[hour:]minute:]second
             f"-L node={self.nodes}",
             f"--mpi proc={self.nproc}",
        ])
        return _call

    def run(self, funcs, single=False, **kwargs):
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
        """
        funcs_fid, kwargs_fid = pickle_function_list(funcs,
                                                     path=self.path.scratch,
                                                     **kwargs)
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
            f"{os.path.join(ROOT_DIR, 'system', 'runscripts', 'run')}",
            f"--funcs {funcs_fid}",
            f"--kwargs {kwargs_fid}",
            f"--environment {self.environs or ''}"
        ])

        # Single-process jobs simply need to replace a few sbatch arguments.
        # Do it AFTER `run_call` has been defined so that subclasses submitting
        # custom run calls can still benefit from this
        if single:
            logger.info("replacing parts of pjsub run call for single "
                        "process job")
            run_call = modify_run_call_single_proc(run_call)

        logger.debug(run_call)

        # Grab the job id (used to monitor job status) from the stdout message
        stdout = subprocess.run(run_call, stdout=subprocess.PIPE,
                                text=True, shell=True).stdout
        job_id = self._stdout_to_job_id(stdout)

        # Monitor the job queue until all jobs have completed, or any one fails
        try:
            status = check_job_status_array(job_id)
        except FileNotFoundError:
            logger.critical(f"cannot access job information through 'sacct', "
                            f"waited 50s with no return, please check job "
                            f"scheduler and log messages")
            sys.exit(-1)

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


def check_job_status_array(job_id):
    """
    Repeatedly check the status of a currently running job using 'sacct'.
    If the job goes into a bad state like 'FAILED', log the failing
    job's id and their states. If all jobs complete nominally, return

    .. note::
        The time.sleep() is critical before querying job status because the 
        system will likely take a second to intitiate jobs so if we 
        `query_job_states` before this has happenend, it will return empty
        lists and cause the function to error out

    :type job_id: str
    :param job_id: main job id to query, returned from the subprocess.run that
    ran the jobs
    :rtype: int
    :return: status of all running jobs. 1 for pass (all jobs COMPLETED). -1 for
        fail (one or more jobs returned failing status)
    :raises FileNotFoundError: if 'sacct' does not return any output for ~1 min.
    """
    logger.info(f"monitoring job status for submitted job: {job_id}")
    while True:
        time.sleep(5)  # give job time to process and also prevent over-query
        job_ids, states = query_job_states(job_id)
        # Sometimes query_job_state() does not return, so we wait again
        if not job_ids or not states:
            continue
        if all([state == "COMPLETED" for state in states]):
            logger.debug("all array jobs returned a 'COMPLETED' state")
            return 1  # Pass
        elif any([check in states for check in BAD_STATES]):  # Any bad states?
            logger.info("atleast 1 system job returned a failing exit code")
            for job_id, state in zip(job_ids, states):
                if state in BAD_STATES:
                    logger.debug(f"{job_id}: {state}")
            return -1  # Fail

def check_job_status_list(job_ids):                                          
    """                                                                          
    Check the status of a list of currently running jobs. This is used for 
    systems that cannot submit array jobs (e.g., Frontera) where we instead
    submit jobs one by one and have to check the status of all those jobs
    together.

    :type job_ids: list of str
    :param job_id: job ID's to query with SACCT. Will be considered one group
        of jobs, who all need to finish successfully otherwise the entire group
        is considered failed
    :rtype: int
    :return: status of all running jobs. 1 for pass (all jobs COMPLETED). -1 for
        fail (one or more jobs returned failing status)
    :raises FileNotFoundError: if 'sacct' does not return any output for ~1 min.
    """                                                                          
    logger.info(f"monitoring job status for {len(job_ids)} submitted jobs")      
                                                                                 
    while True:                                                                  
        time.sleep(10)                                                            
        job_id_list, states = [], []
        for job_id in job_ids:                                                   
            _job_ids, _states = query_job_states(job_id)                         
            job_id_list += _job_ids                                              
            states += _states                                                    
        # Sometimes query_job_state() does not return, so we wait again
        if not states or not job_id_list:
            continue
        if all([state == "COMPLETED" for state in states]):                      
            return 1  # Pass                                                     
        elif any([check in states for check in BAD_STATES]):  # Any bad states?  
            logger.info("atleast 1 system job returned a failing exit code")     
            for job_id, state in zip(job_ids, states):                           
                if state in BAD_STATES:                                          
                    logger.debug(f"{job_id}: {state}")                           
            return -1  # Fail  


def query_job_states(job_id, _recheck=0):
    """
    Queries completion status of an array job by running the SLURM cmd `sacct`
    Available job states are listed here: https://slurm.schedmd.com/sacct.html

    .. note::
        The actual command line call wil look something like this
        $ sacct -nLX -o jobid,state -j 441630
        441630_0    PENDING
        441630_1    COMPLETED

    .. note::
        SACCT flag options are described as follows:
        -L: queries all available clusters, not just the cluster that ran the
            `sacct` call. Used for federated clusters
        -X: supress the .batch and .extern jobnames that are normally returned
            but don't represent that actual running job

    :type job_id: str
    :param job_id: main job id to query, returned from the subprocess.run that
        ran the jobs
    :type rechecks: int
    :param rechecks: Used for recursive calling of the function. It can take 
        time for jobs to be initiated on a system, which may result in the 
        stdout of the 'sacct' command to be empty. In this case we wait and call
        the function again. Rechecks are used to prevent endless loops by 
        putting a stop criteria
    :raises FileNotFoundError: if 'sacct' does not return any output for ~1 min.
    """
    job_ids, job_states = [], []
    cmd = f"sacct -nLX -o jobid,state -j {job_id}"
    stdout = subprocess.run(cmd, stdout=subprocess.PIPE,
                            text=True, shell=True).stdout
    
    # Recursively re-check job state incase the job has not been instantiated 
    # in which cause 'stdout' is an empty string
    if not stdout:
        _recheck += 1
        if _recheck > 10:
            raise FileNotFoundError(f"Cannot access job ID {job_id}")
        time.sleep(10)
        query_job_states(job_id, _recheck)

    # Return the job numbers and respective states for the given job ID
    for job_line in str(stdout).strip().split("\n"):
        if not job_line:
            continue
        job_id, job_state = job_line.split()
        job_ids.append(job_id)
        job_states.append(job_state)

    return job_ids, job_states

def modify_run_call_single_proc(run_call):
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
            run_call = run_call.replace(part, "--array=0-0")
        elif "--ntasks" in part:
            run_call = run_call.replace(part, "--ntasks=1")

    # Append taskid to environment variable, deal with the case where
    # self.par.ENVIRONS is an empty string
    task_id_str = "SEISFLOWS_TASKID=0"
    if not run_call.strip().endswith("--environment"):
        task_id_str = f",{task_id_str}"  # appending to the list of vars

    run_call += task_id_str

    return run_call

