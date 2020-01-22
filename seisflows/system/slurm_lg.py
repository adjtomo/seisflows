#!/usr/bin/env python
"""
This is the subclass seisflows.system.slurm_lg

This class provides the core utilities interaction with HPC systems which run
using Slurm management tools
"""
import os
import math
import sys
import time

from subprocess import check_output
from seisflows.tools import msg, unix
from seisflows.tools.err import ParameterError
from seisflows.tools.tools import call, findpath, saveobj, timestamp
from seisflows.config import custom_import

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class SlurmLg(custom_import('system', 'base')):
    """
    An interface through which to submit workflows, run tasks in serial or
    parallel, and perform other system functions.

    By hiding environment details behind a Python interface layer, these
    classes provide a consistent command set across different computing
    environments.

    Intermediate files are written to a global scratch path PATH.SCRATCH,
    which must be accessible to all compute nodes.

    Optionally, users can provide a local scratch path PATH.LOCAL if each
    compute node has its own local filesystem.

    For important additional information, please see
    http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration
    """
    def check(self):
        """
        Checks parameters and paths
        """
        # Name of job
        if "TITLE" not in PAR:
            setattr(PAR, "TITLE", os.path.basename(os.path.abspath(".")))

        # Time allocated for workflow in minutes
        if "WALLTIME" not in PAR:
            setattr(PAR, "WALLTIME", 30.)

        # Time allocated for each individual task in minutes
        if "TASKTIME" not in PAR:
            setattr(PAR, "TASKTIME", 15.)

        # Number of tasks
        if "NTASK" not in PAR:
            raise ParameterError(PAR, "NTASK")

        # Number of cores per task
        if "NPROC" not in PAR:
            raise ParameterError(PAR, "NPROC")

        # Limit on number of concurrent tasks
        if "NTASKMAX" not in PAR:
            setattr(PAR, "NTASKMAX", 100)

        # Number of cores per node
        if "NODESIZE" not in PAR:
            raise ParameterError(PAR, "NODESIZE")

        # Optional additional SLURM arguments
        if "SLURMARGS" not in PAR:
            setattr(PAR, "SLURMARGS", "")

        # Optional environment variable list VAR1=val1,VAR2=val2,...
        if "ENVIRONS" not in PAR:
            setattr(PAR, "ENVIRONS", "")

        # Level of detail in output messages
        if "VERBOSE" not in PAR:
            setattr(PAR, "VERBOSE", 1)

        # Location where job was submitted
        if "WORKDIR" not in PATH:
            setattr(PATH, "WORKDIR", os.path.abspath("."))

        # Location where output files are written
        if "OUTPUT" not in PATH:
            setattr(PATH, "OUTPUT", os.path.join(PATH.WORKDIR, "output"))

        # Location where temporary files are written
        if "SCRATCH" not in PATH:
            setattr(PATH, "SCRATCH", os.path.join(PATH.WORKDIR, "scratch"))

        # Where system files are written
        if "SYSTEM" not in PATH:
            setattr(PATH, "SYSTEM", os.path.join(PATH.SCRATCH, "system"))

        # Optional local scratch path
        if "LOCAL" not in PATH:
            setattr(PATH, "LOCAL", None)

    def submit(self, workflow):
        """
        Submits workflow as a master job
        """
        # create scratch directories
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # Create output directories
        unix.mkdir(PATH.OUTPUT)
        unix.mkdir(os.path.join(PATH.WORKDIR, "output.slurm"))

        workflow.checkpoint()

        # Submit using sbatch
        submit_call = " ".join([
            f"sbatch {PAR.SLURMARGS}",
            f"--job_name={PAR.TITLE}",
            f"--output={output_log}",
            f"--ntasks-per-node={PAR.NODESIZE}",
            f"--nodes=1",
            f"--time={PAR.WALLTIME:d}",
            os.path.join(findpath("seisflows.system"), "wrappers", "submit"),
            PATH.OUTPUT
        ])

        call(submit_call)

    def run(self, classname, method, *args, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on the
        maui cluster

        Executes classname.method(*args, **kwargs) NTASK times,
        each time on NPROC CPU cores

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        """
        # Checkpoint this individual method before proceeding
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        # Submit job array
        run_call = " ".join([
            "sbatch",
            "{PAR.SLURMARGS}"
            f"--job-name={PAR.TITLE}",
            f"--nodes={math.ceil(PAR.NPROC/float(PAR.NODESIZE)):d}",
            f"--ntasks-per-node={PAR.NODESIZE:d}",
            f"--ntasks={PAR.NPROC:d}",
            f"--time={PAR.TASKTIME:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'output.slurm', '%A_%a')}",
            f"--array=0-{PAR.NTASK-1 & PAR.NTASKMAX}",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}"
        ])

        stdout = check_output(run_call, shell=True)

        # Keep track of job ids
        jobs = self.job_id_list(stdout, PAR.NTASK)

        # Check job array completion status
        while True:
            # Wait a few seconds between queries
            time.sleep(5)
            isdone, jobs = self.job_array_status(classname, method, jobs)
            if isdone:
                return

    def run_single(self, classname, method, *args, **kwargs):
        """
        Runs task a single time

        Executes classname.method(*args, **kwargs) a single time on NPROC
        CPU cores
        """

        # Checkpoint this individual method before proceeding
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        # Submit job array
        run_call = " ".join([
            "sbatch",
            "{PAR.SLURMARGS}"
            f"--job-name={PAR.TITLE}",
            f"--nodes={math.ceil(PAR.NPROC/float(PAR.NODESIZE)):d}",
            f"--ntasks-per-node={PAR.NODESIZE:d}",
            f"--ntasks={PAR.NPROC:d}",
            f"--time={PAR.TASKTIME:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'output.slurm', '%A_%a')}",
            f"--array=0-0",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}"
            f"{SEISFLOWS_TASKID}=0"
        ])

        stdout = check_output(run_call, shell=True)

        # Keep track of job ids
        jobs = self.job_id_list(stdout, PAR.NTASK)

        # Check job array completion status
        while True:
            # Wait a few seconds between queries
            time.sleep(5)
            isdone, jobs = self.job_array_status(classname, method, jobs)
            if isdone:
                return

    def run_ancil(self):
        """
        Special run call to use Python on maui_ancil

        !!! To be implemented by subclass MauiLg !!!
        """
        raise NotImplementedError

    def mpiexec(self):
        """
        Specifies MPI executable used to invoke solver

        :rtype: str
        :return: the MPI exectuabel for a Slurm system
        """
        return 'srun -u '

    def taskid(self):
        """
        Provides a unique identifier for each running task

        !!! TO DO: what type of exception is that? !!!
        :rtype: int
        :return: identifier for a given task
        """
        try:
            return int(os.getenv('SEISFLOWS_TASKID'))
        except:
            return int(os.getenv('SLURM_ARRAY_TASK_ID'))

    def job_array_status(self, classname, method, jobs):
        """
        Determines completion status of job or job array

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        :type jobs: list
        :param jobs: list of jobs currently running
        """
        states = []
        for job in jobs:
            state = self.job_status(job)
            if state in ['TIMEOUT']:
                print(msg.TaskTimeout.format(classname=classname,
                                             method=method, job_id=job,
                                             tasktime=PAR.TASKTIME))
                sys.exit(-1)
            elif state in ['FAILED', 'NODE_FAIL']:
                print(msg.TaskError_SLURM.format(classname=classname,
                                                 method=method, job_id=job))
                sys.exit(-1)
            elif state in ['COMPLETED']:
                states += [1]
            else:
                states += [0]

        isdone = all(states)

        return isdone, jobs

    def job_id_list(self, stdout, ntask):
        """
        Parses job id list from sbatch standard output

        :type stdout: str
        :param stdout: the output of subprocess.check_output()
        :type ntask: int
        :param ntask: number of tasks currently running
        """
        job_id = stdout.split()[-1].strip()
        return [f"{job_id}_{str(ii)}" for ii in range(ntask)]

    def job_status(self, job):
        """
        Queries completion status of a single job

        :type job: str
        :param job: job id to query
        """
        stdout = check_output(
            "sacct -n -o jobid,state -j "+ job.split("_")[0],
            shell=True)

        state = ""
        lines = stdout.strip().split("\n")
        for line in lines:
            if line.split()[0] == job:
                state = line.split()[1]
        return state


