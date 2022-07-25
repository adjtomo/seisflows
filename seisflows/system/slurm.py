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
"""
import os
import sys
import math
import time
import subprocess

from seisflows import logger
from seisflows.system.cluster import Cluster
from seisflows.tools import msg
from seisflows.tools.config import ROOT_DIR


class Slurm(Cluster):
    """
    [system.slurm] Interface for submitting jobs to Simple Linux Utility for
    Resource Management (SLURM) system.

    :type ntask_max: int
    :param ntask_max: limit the number of concurrent tasks in a given array job
    :type slurm_args: str
    :param slurm_args: Any (optional) additional SLURM arguments that will
        be passed to the SBATCH scripts. Should be in the form:
        '--key1=value1 --key2=value2"
    """
    __doc__ = Cluster.__doc__ + __doc__

    def __init__(self, ntask_max=100, slurm_args="",  **kwargs):
        """Slurm-specific setup parameters"""
        super().__init__(**kwargs)

        # Overwrite the existing 'mpiexec'
        if self.mpiexec is None:
            self.mpiexec = "srun -u"
        self.ntask_max = ntask_max
        self.slurm_args = slurm_args

        # Must be overwritten by child class
        self.node_size = None

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        assert(self.node_size is not None), (
            f"Slurm system child classes require defining the `node_size` or "
            f"the number of cores per node inherent to the compute system")

    def submit(self, submit_call=None):
        """
        Submits workflow as a single process master job on a SLURM system

        :type submit_call: str
        :param submit_call: subclasses (e.g., specific SLURM cluster subclasses)
            can overload the sbatch command line input by setting
            submit_call. If set to None, default submit_call will be set here.
        """
        if submit_call is None:
            submit_call = " ".join([
                f"sbatch",
                f"{self.slurm_args or ''}",
                f"--job-name={self.title}",
                f"--output={self.path.output_log}",
                f"--error={self.path.output_log}",
                f"--ntasks-per-node={self.node_size}",
                f"--nodes=1",
                f"--time={self.walltime:d}",
                f"{os.path.join(ROOT_DIR, 'system', 'runscripts', 'submit')}",
                f"--output {self.path.output}"
            ])

        logger.debug(submit_call)
        super().submit(submit_call=submit_call)

    def run(self, funcs, single=False, run_call=None, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on a SLURM
        cluster. Executes classname.method(*args, **kwargs) `NTASK` times,
        each time on `NPROC` CPU cores

        :type funcs: list of methods
        :param funcs: a list of functions that should be run in order. All
            kwargs passed to run() will be passed into the functions.
        :type single: bool
        :param single: run a single-process, non-parallel task, such as
            smoothing the gradient, which only needs to be run by once.
            This will change how the job array and the number of tasks is
            defined, such that the job is submitted as a single-core job to
            the system.
        :type run_call: str
        :param run_call: subclasses (e.g., specific SLURM cluster subclasses)
            can overload the sbatch command line input by setting
            run_call. If set to None, default run_call will be set here.
        """
        funcs_fid, kwargs_fid = self._pickle_func_list(funcs, **kwargs)
        logger.info(f"running functions {[_.__name__ for _ in funcs]} on "
                    f"system {self.ntask} times")

        # Default sbatch command line input, can be overloaded by subclasses
        # Copy-paste this default run_call and adjust accordingly for subclass
        if run_call is None:
            run_call = " ".join([
                "sbatch",
                f"{self.slurm_args or ''}",
                f"--job-name={self.title}",
                f"--nodes={math.ceil(self.nproc/float(self.node_size)):d}",
                f"--ntasks-per-node={self.node_size:d}",
                f"--ntasks={self.nproc:d}",
                f"--time={self.tasktime:d}",
                f"--output={os.path.join(self.path.log_files, '%A_%a')}",
                f"--array=0-{self.ntask-1}%{self.ntask_max}",
                f"--parsable",  # keeps stdout cleaner
                f"{os.path.join(ROOT_DIR, 'system', 'runscripts', 'run')}",
                f"--funcs {funcs_fid}",
                f"--kwargs {kwargs_fid}",
                f"--environment {self.environs or ''}"
            ])
        logger.debug(run_call)

        # Single-process jobs simply need to replace a few sbatch arguments.
        # Do it AFTER `run_call` has been defined so that subclasses submitting
        # custom run calls can still benefit from this
        if single:
            logger.info("replacing parts of sbatch run call for single "
                        "process job")
            run_call = _modify_run_call_single_proc(run_call)

        # Stdout will be job number (e.g., 1234). Federated clusters will return
        # job # and cluster name (e.g., 1234;Cluster1). We only want job #
        job_id = subprocess.run(run_call, stdout=subprocess.PIPE,
                                text=True, shell=True).stdout
        job_id = str(job_id).split(";")[0]

        # Monitor the job queue until all jobs have completed, or any one fails
        status = check_job_status(job_id)
        if status == -1:  # Failed job
            logger.critical(
                msg.cli(f"Stopping workflow. Please check logs for details.",
                        items=[f"TASKS:   {[_.__name_ for _ in funcs]}",
                               f"SBATCH:  {run_call}"],
                        header="slurm run error", border="=")
            )
            sys.exit(-1)
        else:
            logger.info(f"tasks finished successfully")


def check_job_status(job_id):
    """
    Repeatedly check the status of a currently running job using 'sacct'.
    If the job goes into a bad state like 'FAILED', log the failing
    job's id and their states. If all jobs complete nominally, return

    :type job_id: str
    :param job_id: main job id to query, returned from the subprocess.run that
    ran the jobs
    :rtype: int
    :return: status of all running jobs. 1 for pass (all jobs COMPLETED). -1 for
        fail (one or more jobs returned failing status)
    """
    bad_states = ["TIMEOUT", "FAILED", "NODE_FAIL",
                  "OUT_OF_MEMORY", "CANCELLED"]
    while True:
        job_ids, states = query_job_states(job_id)
        if [state == "COMPLETED" for state in states]:
            return 1  # Pass
        elif any([check in states for check in bad_states]):  # Any bad states?
            logger.info("atleast 1 system job returned a failing exit code")
            for job_id, state in zip(job_ids, states):
                if state in bad_states:
                    logger.debug(f"{job_id}: {state}")
            return -1  # Fail
        else:
            time.sleep(5)  # Don't query 'sacct' command too often


def query_job_states(job_id):
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
    """
    job_ids, job_states = [], []
    cmd = f"sacct -nLX -o jobid,state -j {job_id}"
    stdout = subprocess.run(cmd, stdout=subprocess.PIPE,
                            text=True, shell=True).stdout
    for job_line in str(stdout).strip().split("\n"):
        job_id, job_state = job_line.split()
        job_ids.append(job_id)
        job_states.append(job_state)

    return job_ids, job_states


def _modify_run_call_single_proc(run_call):
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
            run_call.replace(part, "--array=0-0")
        elif "--ntasks" in part:
            run_call.replace(part, "--ntasks=1")

    # Append taskid to environment variable, deal with the case where
    # self.par.ENVIRONS is an empty string
    task_id_str = "SEISFLOWS_TASKID=0"
    if not run_call.strip().endswith("--environment"):
        task_id_str = f",{task_id_str}"  # appending to the list of vars

    run_call += task_id_str

    return run_call

