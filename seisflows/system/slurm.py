#!/usr/bin/env python3
"""
The Simmple Linux Utility for Resource Management (SLURM) is a commonly used
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
import logging
import subprocess

from seisflows.tools import msg
from seisflows.config import ROOT_DIR, custom_import, SeisFlowsPathsParameters

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Slurm(custom_import("system", "cluster")):
    """
    Generalized interface for submitting jobs to and interfacing with a SLURM
    workload management system.
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        super().__init__()

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        sf.par("MPIEXEC", required=False, default="srun -u", par_type=str,
               docstr="Function used to invoke executables on the system. "
                      "For example 'srun' on SLURM systems, or './' on a "
                      "workstation. If left blank, will guess based on the "
                      "system.")

        sf.par("NTASKMAX", required=False, default=100, par_type=int,
               docstr="Limit on the number of concurrent tasks in array")

        sf.par("NODESIZE", required=True, par_type=int,
               docstr="The number of cores per node defined by the system")

        sf.par("SLURMARGS", required=False, default="", par_type=str,
               docstr="Any optional, additional SLURM arguments that will be "
                      "passed to the SBATCH scripts")

        return sf

    def submit(self, submit_call=None):
        """
        Submits workflow as a single process master job

        :type workflow: module
        :param workflow:
        :type submit_call: str
        :param submit_call: subclasses (e.g., specific SLURM cluster subclasses)
            can overload the sbatch command line input by setting
            submit_call. If set to None, default submit_call will be set here.
        """
        if submit_call is None:
            submit_call = " ".join([
                f"sbatch",
                f"{PAR.SLURMARGS or ''}",
                f"--job-name={PAR.TITLE}",
                f"--output={self.output_log}",
                f"--error={self.error_log}",
                f"--ntasks-per-node={PAR.NODESIZE}",
                f"--nodes=1",
                f"--time={PAR.WALLTIME:d}",
                f"{os.path.join(ROOT_DIR, 'scripts', 'submit')}",
                "--output {PATH.OUTPUT}"
            ])
            self.logger.debug(submit_call)

        super().submit(submit_call)

    def run(self, classname, method, single=False, run_call=None, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on a SLURM
        cluster. Executes classname.method(*args, **kwargs) `NTASK` times,
        each time on `NPROC` CPU cores

        .. note::
            The actual CLI call structure looks something like this
            $ sbatch --args scripts/run OUTPUT class method environs

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
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
        self.checkpoint(PATH.OUTPUT, classname, method, kwargs)

        # Default sbatch command line input, can be overloaded by subclasses
        # Copy-paste this default run_call and adjust accordingly for subclass
        if run_call is None:
            run_call = " ".join([
                "sbatch",
                f"{PAR.SLURMARGS or ''}",
                f"--job-name={PAR.TITLE}",
                f"--nodes={math.ceil(PAR.NPROC/float(PAR.NODESIZE)):d}",
                f"--ntasks-per-node={PAR.NODESIZE:d}",
                f"--ntasks={PAR.NPROC:d}",
                f"--time={PAR.TASKTIME:d}",
                f"--output={os.path.join(PATH.WORKDIR, 'logs', '%A_%a')}",
                f"--array=0-{PAR.NTASK-1 % PAR.NTASKMAX}",
                f"{os.path.join(ROOT_DIR, 'scripts', 'run')}",
                f"--output {PATH.OUTPUT}",
                f"--classname {classname}",
                f"--funcname {method}",
                f"--environment {PAR.ENVIRONS or ''}"
            ])
            self.logger.debug(run_call)

        # Single-process jobs simply need to replace a few sbatch arguments.
        # Do it AFTER `run_call` has been defined so that subclasses submitting
        # custom run calls can still benefit from this
        if single:
            self.logger.info("replacing parts of sbatch run call for single "
                             "process job")
            run_call = _modify_run_call_single_proc(run_call)

        # The standard response from SLURM when submitting jobs
        # is something like 'Submitted batch job 441636', we want job number
        stdout = subprocess.run(run_call, stdout=subprocess.PIPE,
                                text=True, shell=True).stdout

        # Continuously check for job completion on ALL running array jobs
        job_ids = job_id_list(stdout, single)
        is_done = False
        count = 0
        bad_states = ["TIMEOUT", "FAILED", "NODE_FAIL", "OUT_OF_MEMORY",
                      "CANCELLED"]
        while not is_done:
            # Wait a bit to avoid rapidly querying sacct
            time.sleep(5)
            is_done, states = job_array_status(job_ids)
            # EXIT CONDITION: if any of the jobs provide job failure codes
            if not is_done:
                for i, state in enumerate(states):
                    # Sometimes states can be something like 'CANCELLED+', so
                    # we can't do exact string matching, check partial matches
                    if any([check in state for check in bad_states]):
                        print(msg.cli((f"Stopping workflow for {state} job. "
                                       f"Please check log file for details."),
                                      items=[f"TASK:    {classname}.{method}",
                                             f"TASK ID: {job_ids[i]}",
                                             f"LOG:     logs/{job_ids[i]}",
                                             f"SBATCH:  {run_call}"],
                                      header="slurm run error", border="="))
                        sys.exit(-1)
                        # WAIT CONDITION: if sacct is not working, we'll get stuck in a loop
            if "UNDEFINED" in states:
                count += 1
                # Every 10 counts, warn the user this is unexpected behavior
                if not count % 10:
                    job_id = job_ids[states.index("UNDEFINED")]
                    self.logger.warning(f"SLURM command 'sacct {job_id}' has "
                                        f"returned unexpected response {count} "
                                        f"times. This job may have failed "
                                        f"unexpectedly. Consider checking "
                                        f"manually")

        self.logger.info(f"Task {classname}.{method} finished successfully")

    def taskid(self):
        """
        Provides a unique identifier for each running task

        :rtype: int
        :return: identifier for a given task
        """
        # If not set, this environment variable will return None
        sftaskid = os.getenv("SEISFLOWS_TASKID")

        if sftaskid is None:
            sftaskid = os.getenv("SLURM_ARRAY_TASK_ID")
            if sftaskid is None:
                print(msg.cli("system.taskid() environment variable not found. "
                              "Assuming DEBUG mode and returning taskid==0. "
                              "If not DEBUG mode, please check SYSTEM.run()",
                              header="warning", border="="))
                sftaskid = 0

        return int(sftaskid)


def job_id_list(stdout, single):
    """
    Parses job id list from sbatch standard output. Stdout typically looks
    like: 'Submitted batch job 441636', but if submitting jobs cross-cluster
    (e.g., like on Maui), stdout might be:
    'Submitted batch job 441636 on cluster Maui'

    .. note::
        In order to find the job number, we just scan each word in stdout
        until we find the number, ASSUMING that there is only one number in
        the string

    TODO Should failing to return job_id break in reasonable way?

    The output job arrays will look something like:
    [44163_0, 44163_1, ..., 44163_PAR.NTASK]

    :type stdout: str
    :param stdout: the text response from running 'sbatch' on SLURM, which
        should be returned by subprocess.run(stdout=PIPE)
    :type single: bool
    :param single: if running a single process job, returns a list of length
        1 with a single job id, else returns a list of length PAR.NTASK
        for all arrayed jobs
    :rtype: list
    :return: a list of array jobs that should be currently running
    """
    if single:
        ntask = 1
    else:
        ntask = PAR.NTASK

    # Splitting e.g.,: 'Submitted batch job 441636\n'
    for part in stdout.strip().split():
        try:
            # The int will keep throwing ValueError until we find the num
            job_id = int(part)
            break
        except ValueError:
            continue
    return [f"{job_id}_{i}" for i in range(ntask)]


def job_array_status(job_ids):
    """
    Determines current status of job or job array

    :type job_ids: list
    :param job_ids: list of SLURM job id numbers to check completion of
        Will not return unless all jobs have completed
    :rtype is_done: bool
    :return is_done: True if all jobs in the array have been completed
    :rtype states: list
    :return states: list of states returned from sacct
    """
    states = []
    for job_id in job_ids:
        state = check_job_state(job_id)
        states.append(state.upper())

    # All array jobs must be completed to return is_done == True
    is_done = all([state.upper() == "COMPLETED" for state in states])

    return is_done, states


def check_job_state(job_id):
    """
    Queries completion status of a single job by running:
        $ sacct -nL -o jobid,state -j {job_id}

        # Example outputs from this sacct command
        # JOB_ID    STATUS
        441630_0  PENDING  # array job will have the array number
        441630    COMPLETED  # if --array=0-0, jobs will not have suffix
        441628.batch    COMPLETED  # we don't want to check these

    Available job states: https://slurm.schedmd.com/sacct.html

    .. note::
        -L flag in sacct queries all available clusters, not just the
        cluster that ran the `sacct` call
        -X supress the .batch and .extern jobname

    :type job: str
    :param job: job id to query
    """
    cmd = f"sacct -nLX -o jobid,state -j {job_id}"
    stdout = subprocess.run(cmd, stdout=subprocess.PIPE,
                            text=True, shell=True).stdout

    # Undefined status will be retured if we cannot match the job id with
    # the sacct output
    state = "UNDEFINED"
    lines = stdout.strip().split("\n")
    for line in lines:
        # expecting e.g., 441628  COMPLETED
        try:
            job_id_check, state = line.split()
        # str.split() will throw ValueError on non-matching strings
        except ValueError:
            continue
        # Use in to allow for array jobs to match job ids
        if job_id in job_id_check:
            break

    return state


def _modify_run_call_single_proc(run_call):
    """
    Modifies a SLURM SBATCH command to use only 1 processor as a single run

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
    # PAR.ENVIRONS is an empty string
    task_id_str = "SEISFLOWS_TASKID=0"
    if not run_call.strip().endswith("--environment"):
        task_id_str = f",{task_id_str}"  # appending to the list of vars

    run_call += task_id_str

    return run_call

