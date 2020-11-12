#!/usr/bin/env python
"""
This is the subclass seisflows.system.maui_lg
This class provides the core utilities interaction with HPC systems which must
be overloaded by subclasses
"""
import os
import sys
import math
import time
from subprocess import check_output, call, CalledProcessError
from seisflows.tools.tools import call, findpath
from seisflows.config import custom_import, SeisFlowsPathsParameters

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class Maui(custom_import("system", "slurm_lg")):
    """
    System interface for the New Zealand Tomography problem

    Inversions are run on New Zealand eScience Infrascructure (NeSI) HPCs
    Heavy simulation work is run on Maui
    Preprocessing tasks, via Pyatoa, are run on Maui_ancil, the ancillary
    cluster attached to Maui.

    Both clusters are run with the Slurm system, and
    so MauiLG inherits attributes from `slurm_lg` system
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        sf.par("ACCOUNT", required=True, par_type=str,
               docstr="The account name to submit jobs under")

        sf.par("NODES", required=True, par_type=int,
               docstr="The number of nodes to use per job. Rough estimate "
                      "would be NPROC//NODESIZE")

        sf.par("CPUS_PER_TASK", required=False, default=1, par_type=int,
               docstr="Multiple CPUS per task allows for multithreading jobs")

        sf.par("MAIN_CLUSTER", required=False, default="maui", par_type=str,
               docstr="Name of main cluster for job submission")

        sf.par("MAIN_PARTITION", required=False, default="nesi_research",
               par_type=str, docstr="Name of partition on main cluster")

        sf.par("ANCIL_CLUSTER", required=False, default="maui_ancil",
               par_type=str,
               docstr="Name of ancillary cluster for prepost tasks")

        sf.par("ANCIL_PARTITION", required=False, default="nesi_prepost",
               par_type=str,
               docstr="Name of ancillary partition for prepost tasks")

        sf.par("ANCIL_TASKTIME", required=False, default="null", par_type=float,
               docstr="Tasktime for prepost jobs on ancillary nodes")

        sf.par("NODESIZE", required=False, default=40, par_type=int,
               docstr="The number of cores per node defined by the system")

        sf.par("MPIEXEC", required=False, default="srun", par_type=str,
               docstr="Function used to invoke parallel executables")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

        assert(PAR.NODESIZE == 40), f"Maui NODESIZE is physically set to 40"

    def submit(self, workflow):
        """
        Overwrites seisflows.workflow.slurm_lg.submit()

        Submits master job workflow to maui_ancil cluster with a more in-depth
        logging system that saves old output logs and names logs based on
        the job name

        Note:
            The master job must be run on maui_ancil because Maui does
            not have the ability to run the command "sacct"
        """
        output_log, error_log = self.setup()
        workflow.checkpoint()

        # Submit to maui_ancil
        submit_call = " ".join([
            f"sbatch {PAR.SLURMARGS}",
            f"--account={PAR.ACCOUNT}",
            f"--cluster={PAR.ANCIL_CLUSTER}",
            f"--partition={PAR.ANCIL_PARTITION}",
            f"--job-name=main_{PAR.TITLE}",  # main_ prefix means master
            f"--output={output_log}-%A.log",
            f"--error={error_log}-%A.log",
            f"--ntasks=1",
            f"--cpus-per-task=1",
            f"--time={PAR.WALLTIME:d}",
            os.path.join(findpath("seisflows.system"), "wrappers", "submit"),
            PATH.OUTPUT
        ])
        call(submit_call)

    def run(self, classname, method, scale_tasktime=1, *args, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on the
        maui cluster

        Executes classname.method(*args, **kwargs) NTASK times,
        each time on NPROC CPU cores

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        :type scale_tasktime: int
        :param scale_tasktime: a way to get over the hard-set tasktime, because
            some tasks take longer (e.g. smoothing), but you don't want these
            to set the tasktimes for all other tasks. This lets you scale the
            time of specific tasks by PAR.TASKTIME * scale_tasktime
        """
        # Checkpoint this individual method before proceeding
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        run_call = " ".join([
            "sbatch",
            f"{PAR.SLURMARGS}",
            f"--account={PAR.ACCOUNT}",
            f"--job-name={PAR.TITLE}",
            f"--clusters={PAR.MAIN_CLUSTER}",
            f"--partition={PAR.MAIN_PARTITION}",
            f"--cpus-per-task={PAR.CPUS_PER_TASK}",
            f"--nodes={PAR.NODES:d}",
            f"--ntasks={PAR.NPROC:d}",
            f"--time={PAR.TASKTIME * scale_tasktime:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'output.logs', '%A_%a')}",
            f"--array=0-{PAR.NTASK-1 % PAR.NTASKMAX}",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}"
        ])

        stdout = check_output(run_call, shell=True)
        
        # Keep track of Job IDS
        jobs = self.job_id_list(stdout, PAR.NTASK)

        # Check for job completion status
        check_status_error = 0
        while True:
            # Wait a few seconds between queries
            time.sleep(5)

            # Occassionally connections using 'sacct' are refused leading to job
            # failure. Wrap in a try-except and allow a handful of failures 
            # incase the failure was a one-off connection problem
            try:
                isdone, jobs = self.job_array_status(classname, method, jobs)
            except CalledProcessError:
                check_status_error += 1
                if check_status_error >= 10:
                    print("check job status with sacct failed 10 times")
                    sys.exit(-1)
                pass
            if isdone:
                return
        
    def run_single(self, classname, method, scale_tasktime=1, *args, **kwargs):
        """
        Runs task a single time

        Executes classname.method(*args, **kwargs) a single time on NPROC
        CPU cores

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        :type scale_tasktime: int
        :param scale_tasktime: a way to get over the hard-set tasktime, because
            some tasks take longer (e.g. smoothing), but you don't want these
            to set the tasktimes for all other tasks. This lets you scale the
            time of specific tasks by PAR.TASKTIME * scale_tasktime
        """
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        # Submit job
        run_call = " ".join([
            "sbatch",
            f"{PAR.SLURMARGS}",
            f"--account={PAR.ACCOUNT}",
            f"--job-name={PAR.TITLE}",
            f"--clusters={PAR.MAIN_CLUSTER}",
            f"--partition={PAR.MAIN_PARTITION}",
            f"--cpus-per-task={PAR.CPUS_PER_TASK}",
            f"--nodes={PAR.NODES:d}",
            f"--ntasks={PAR.NPROC:d}",
            f"--time={PAR.TASKTIME * scale_tasktime:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'output.logs', '%A_%a')}",
            f"--array=0-0",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}",
            f"SEISFLOWS_TASKID=0"
        ])
        
        stdout = check_output(run_call, shell=True)
        
        # Keep track of job ids
        jobs = self.job_id_list(stdout, 1)

        # Check for job completion status
        check_status_error = 0
        while True:
            # Wait a few seconds between queries
            time.sleep(5)

            # Occassionally connections using 'sacct' are refused leading to job
            # failure. Wrap in a try-except and allow a handful of failures
            # incase the failure was a one-off connection problem
            try:
                isdone, jobs = self.job_array_status(classname, method, jobs)
            except CalledProcessError:
                check_status_error += 1
                if check_status_error >= 10:
                    print("check job status with sacct failed 10 times")
                    sys.exit(-1)
                pass
            if isdone:
                return

    def run_ancil(self, classname, method, *args, **kwargs):
        """
        Runs task a single time. For Maui this is run on maui ancil
        and also includes some extra arguments for eval_func

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        """
        # Set the tasktime required for ancillary tasks
        if PAR.ANCIL_TASKTIME is None:
            ANCIL_TASKTIME = PAR.TASKTIME
        else:
            ANCIL_TASKTIME = PAR.ANCIL_TASKTIME

        # Checkpoint this individual method before proceeding
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        run_call = " ".join([
            "sbatch",
            f"{PAR.SLURMARGS}",
            f"--account={PAR.ACCOUNT}",
            f"--job-name={PAR.TITLE}",
            f"--clusters={PAR.ANCIL_CLUSTER}",
            f"--partition={PAR.ANCIL_PARTITION}",
            f"--cpus-per-task={PAR.CPUS_PER_TASK}",
            f"--time={ANCIL_TASKTIME:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'output.logs', '%A_%a')}",
            f"--array=0-{PAR.NTASK-1 % PAR.NTASKMAX}",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}",
        ])

        stdout = check_output(run_call, shell=True)

        # Keep track of job ids
        jobs = self.job_id_list(stdout, PAR.NTASK)

        # Check for job completion status
        check_status_error = 0
        while True:
            # Wait a few seconds between queries
            time.sleep(5)

            # Occassionally connections using 'sacct' are refused leading to job
            # failure. Wrap in a try-except and allow a handful of failures
            # incase the failure was a one-off connection problem
            try:
                isdone, jobs = self.job_array_status(classname, method, jobs)
            except CalledProcessError:
                check_status_error += 1
                if check_status_error >= 10:
                    print("check job status with sacct failed 10 times")
                    sys.exit(-1)
                pass
            if isdone:
                return

    def job_id_list(self, stdout, ntask):
        """
        Overwrite seisflows.system.workflow.slurm_log.job_id_list()

        Parses job id list from sbatch standard output. 
        Decode class bytes to str using UTF-8 from subprocess.check_output()

        Note:
            Submitting jobs across clusters on Maui means the phrase
            "on cluster X" gets appended to the end of stdout and the job is no
            longer stdout().split()[-1]. Instead, scan through stdout and try
            to find the number using float() to break on words.

        :type stdout: str or bytes
        :param stdout: the output of subprocess.check_output()
        :type ntask: int
        :param ntask: number of tasks currently running
        """    
        if isinstance(stdout, bytes):
            stdout = stdout.decode("UTF-8")

        for parts in stdout.split():
            try:
                job_id = parts.strip()
                _ = float(job_id)  # this can raise a ValueError
                return [f"{job_id}_{str(ii)}" for ii in range(ntask)]
            except ValueError:
                continue


    

