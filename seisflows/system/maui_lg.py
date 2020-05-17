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
from glob import glob
from subprocess import check_output, call, CalledProcessError

from seisflows.tools import unix
from seisflows.tools.tools import call, findpath
from seisflows.config import custom_import

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class MauiLg(custom_import('system', 'slurm_lg')):
    """
    System interface for the New Zealand Tomography problem

    Inversions are run on New Zealand eScience Infrascructure (NeSI) HPCs
    Heavy simulation work is run on Maui
    Preprocessing tasks, via Pyatoa, are run on Maui_ancil, the ancillary
    cluster attached to Maui.

    Both clusters are run with the Slurm system, and
    so MauiLG inherits attributes from `slurm_lg` system
    """
    def check(self):
        """
        Checks parameters and paths
        """
        # Run SlurmLG checks first
        super(MauiLg, self).check()

        # NeSI Nodesize is hard set to 40
        if PAR.NODESIZE != 40:
            print("Maui must have a nodesize of 40, overwriting user set")
            setattr(PAR, "NODESIZE", 40)

        # How to invoke executables
        if "MPIEXEC" not in PAR:
            setattr(PAR, "MPIEXEC", "srun")
        
        # NeSI cluster Maui and Maui Anciliary Node specific variables
        if "ACCOUNT" not in PAR:
            raise Exception("Must specify the 'ACCOUNT' to submit jobs under")

        # Main cluster name
        if "MAIN_CLUSTER" not in PAR:
            setattr(PAR, "MAIN_CLUSTER", "maui")

        # Specific partition of the main cluster
        if "MAIN_PARTITION" not in PAR:
            setattr(PAR, "MAIN_PARTITION", "nesi_research")

        # Ancilary cluster name for pre-processing
        if "ANCIL_CLUSTER" not in PAR:
            setattr(PAR, "ANCIL_CLUSTER", "maui_ancil")

        # Specific partition of ancilary cluster
        if "ANCIL_PARTITION" not in PAR:
            setattr(PAR, "ANCIL_PARTITION", "nesi_prepost")

        # If preprocessing tasks are much less than PAR.TASKTIME, it can be
        # useful to manually set a shorter tasktime
        if "ANCIL_TASKTIME" not in PAR:
            setattr(PAR, "ANCIL_TASKTIME", PAR.TASKTIME)

        # If number of nodes not given, automatically calculate.
        # if the "nomultithread" hint is given, the number of nodes will need 
        # to be manually set
        if "NODES" not in PAR:
            setattr(PAR, "NODES", math.ceil(PAR.NPROC/float(PAR.NODESIZE)))

        # Allows for multithreading
        if "CPUS_PER_TASK" not in PAR:
            setattr(PAR, "CPUS_PER_TASK", 1)

        # Specfem3D has deprecated OpenMP but the option to submit with it stays
        if "WITH_OPENMP" not in PAR:
            setattr(PAR, "WITH_OPENMP", False)

        # Make sure that y * z = c * x; where x=nodes, y=ntasks= z=cpus-per-task
        # and c=number of cores per node, if using OpenMP
        if PAR.WITH_OPENMP:
            assert(PAR.NPROC * PAR.CPUS_PER_TASK == PAR.NODESIZE * PAR.NODES)

    def submit(self, workflow):
        """
        Overwrites seisflows.workflow.maui_lg.submit()

        Submits master job workflow to maui_ancil cluster

        Note:
            The master job must be run on maui_ancil because Maui does
            not have the ability to run the command "sacct"
        """
        # Create scratch directories
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # Create output directories
        unix.mkdir(PATH.OUTPUT)
        unix.mkdir(os.path.join(PATH.WORKDIR, "output.slurm"))
        unix.mkdir(os.path.join(PATH.WORKDIR, "logs"))

        # If a scratch directory is made outside the working directory
        if not os.path.exists('./scratch'):
            unix.ln(PATH.SCRATCH, os.path.join(PATH.WORKDIR, "scratch"))
        
        output_log = os.path.join(PATH.WORKDIR, "output")
        error_log = os.path.join(PATH.WORKDIR, "error")
    
        # If resuming, move old log files to keep them out of the way
        for log in [output_log, error_log]:
            unix.mv(src=glob(os.path.join(f"{log}*.log")), 
                    dst=os.path.join(PATH.WORKDIR, "logs")
                    )
        
        # Copy the parameter.yaml file into the log directoroy
        par_copy = f"parameters_{PAR.BEGIN}-{PAR.END}.yaml"
        unix.cp(src="parameters.yaml", 
                dst=os.path.join(PATH.WORKDIR, "logs", par_copy)
                )
                    
        workflow.checkpoint()

        # Submit to maui_ancil
        submit_call = " ".join([
            f"sbatch {PAR.SLURMARGS}",
            f"--account={PAR.ACCOUNT}",
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

    def prep_openmp(self, subprocess_call):
        """
        OpenMP requires some initial exports before sbatch command can be run
        This function will append these to the 'check_output' call that 'run'
        and 'run_single' use

        :type subprocess_call: str
        :param subprocess_call: the string that is passed to check_output
        :rtype: str
        :return: a prepended call with the correct export statements
        """
        prepended_call = " ".join([
                    "export OMP_NUM_THREADS={PAR.CPUS_PER_TASK};",
                    "export OMP_PROC_BIND=true;",
                    "export OMP_PLACES=cores;",
                    subprocess_call
                    ])
        
        return prepended_call

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
            f"--output={os.path.join(PATH.WORKDIR, 'output.slurm', '%A_%a')}",
            f"--array=0-{PAR.NTASK-1 % PAR.NTASKMAX}",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}"
        ])

        if PAR.WITH_OPENMP:
            run_call = self.prep_openmp(run_call) 

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
            f"--output={os.path.join(PATH.WORKDIR, 'output.slurm', '%A_%a')}",
            f"--array=0-0",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}",
            f"SEISFLOWS_TASKID=0"
        ])

        if PAR.WITH_OPENMP:
            run_call = self.prep_openmp(run_call)
        
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
            f"--time={PAR.ANCIL_TASKTIME:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'output.slurm', '%A_%a')}",
            f"--array=0-{PAR.NTASK-1 % PAR.NTASKMAX}",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}",
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

    def job_status(self, job):
        """
        Overwrite seisflows.system.workflow.slurm_log.job_status()

        Queries completion status of a single job

        The added -L flag to `sacct` to query all clusters

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

    

