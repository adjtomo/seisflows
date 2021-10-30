#!/usr/bin/env python
"""
This is the subclass seisflows.system.slurm_lg

This class provides the core utilities interaction with HPC systems which run
using Slurm management tools

..note::
    The main development system for SeisFlows3 used SLURM. Therefore the other
    system supers will not be up to date until access to those systems are
    granted. This rosetta stone, for converting from SLURM to other workload
    management tools will be useful: https://slurm.schedmd.com/rosetta.pdf
"""
import os
import math
import sys
import time

from warnings import warn
from subprocess import check_output
from seisflows3.tools import msg, unix
from seisflows3.tools.err import ParameterError
from seisflows3.tools.tools import call, findpath
from seisflows3.config import custom_import, SeisFlowsPathsParameters

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class SlurmLg(custom_import("system", "base")):
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
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        # Define the Parameters required by this module
        sf.par("NTASKMAX", required=False, default=100, par_type=int,
               docstr="Limit on the number of concurrent tasks in array")

        sf.par("NODESIZE", required=True, par_type=int,
               docstr="The number of cores per node defined by the system")

        sf.par("SLURMARGS", required=False, default="", par_type=str,
               docstr="Any optional, additional SLURM arguments that will be "
                      "passed to the SBATCH scripts")

        sf.par("ENVIRONS", required=False, default="", par_type=str,
               docstr="Optional environment variables to be provided in the"
                      "following format VAR1=var1,VAR2=var2...")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

    def submit(self, workflow):
        """
        Submits workflow as a master job
        """
        output_log, error_log = self.setup()
        workflow.checkpoint()

        # Submit using sbatch
        submit_call = " ".join([
            f"sbatch", 
            f"{PAR.SLURMARGS or ''}",
            f"--job-name={PAR.TITLE}",
            f"--output={output_log}-%A.log",
            f"--error={error_log}-%A.log",
            f"--ntasks-per-node={PAR.NODESIZE}",
            f"--nodes=1",
            f"--time={PAR.WALLTIME:d}",
            os.path.join(findpath("seisflows3.system"), "wrappers", "submit"),
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
            f"{PAR.SLURMARGS or ''}",
            f"--job-name={PAR.TITLE}",
            f"--nodes={math.ceil(PAR.NPROC/float(PAR.NODESIZE)):d}",
            f"--ntasks-per-node={PAR.NODESIZE:d}",
            f"--ntasks={PAR.NPROC:d}",
            f"--time={PAR.TASKTIME:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'output.logs', '%A_%a')}",
            f"--array=0-{PAR.NTASK-1 % PAR.NTASKMAX}",
            f"{os.path.join(findpath('seisflows3.system'), 'wrappers', 'run')}",
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
            f"{PAR.SLURMARGS or ''}",
            f"--job-name={PAR.TITLE}",
            f"--nodes={math.ceil(PAR.NPROC/float(PAR.NODESIZE)):d}",
            f"--ntasks-per-node={PAR.NODESIZE:d}",
            f"--ntasks={PAR.NPROC:d}",
            f"--time={PAR.TASKTIME:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'output.logs', '%A_%a')}",
            f"--array=0-0",
            f"{os.path.join(findpath('seisflows3.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}"
            f"SEISFLOWS_TASKID=0"
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

        :rtype: int
        :return: identifier for a given task
        """
        try:
            return int(os.getenv('SEISFLOWS_TASKID'))
        except TypeError:
            try:
                # TypeError thrown when 'SEISFLOWS_TASKID' is returned as None
                return int(os.getenv('SLURM_ARRAY_TASK_ID'))    
            except TypeError:
                # This should only be thrown during debug processes when no 
                # job have been submitted but 'taskid' may be called from the
                # debugger
                warn("No task ID; returning index 0. "
                     "Assuming debug environment, if not please check system.")
                return 0
                
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
        Decode class bytes to str using UTF-8

        :type stdout: str
        :param stdout: the output of subprocess.check_output()
        :type ntask: int
        :param ntask: number of tasks currently running
        """
        if isinstance(stdout, bytes):
            stdout = stdout.decode("UTF-8")

        job_id = stdout.split()[-1].strip()
        return [f"{job_id}_{str(ii)}" for ii in range(ntask)]

    def job_status(self, job):
        """
        Queries completion status of a single job

            -L flag in sacct queries all available clusters, not just the
            cluster that ran the `sacct` call

        :type job: str
        :param job: job id to query
        """
        stdout = check_output(
            "sacct -nL -o jobid,state -j " + job.split("_")[0],
            shell=True)

        if isinstance(stdout, bytes):
            stdout = stdout.decode("UTF-8")

        state = ""
        lines = stdout.strip().split("\n")
        for line in lines:
            if line.split()[0] == job:
                state = line.split()[1]
        return state


