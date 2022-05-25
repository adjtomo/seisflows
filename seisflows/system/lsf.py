#!/usr/bin/env python3
"""
This is the subclass seisflows.system.lsf_lg

This class provides the core utilities interaction with HPC systems which run
using the Platform Load Sharing Facility (LSF) workload management platform.
"""
import os
import sys
import time
import logging
import subprocess

from seisflows.tools import msg, unix
from seisflows.tools.wrappers import findpath
from seisflows.config import custom_import, SeisFlowsPathsParameters


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class Lsf(custom_import("system", "cluster")):
    """
    An interface through which to submit workflows, run tasks in serial or
    parallel, and perform other system functions.

    By hiding environment details behind a python interface layer, these
    classes provide a consistent command set across different computing
    environments.

    Intermediate files are written to a global scratch path PATH.SCRATCH,
    which must be accessible to all compute nodes.

    Optionally, users can provide a local scratch path PATH.LOCAL if each
    compute node has its own local filesystem.

    For important additional information, please see
    http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration
    """
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
        Checks parameters and paths
        """
        sf = SeisFlowsPathsParameters(super().required)

        sf.par("MPIEXEC", required=False, default="mpiexec", par_type=str,
               docstr="Function used to invoke executables on the system. "
                      "For example 'srun' on SLURM systems, or './' on a "
                      "workstation. If left blank, will guess based on the "
                      "system.")

        # Define the Parameters required by this module
        sf.par("NTASKMAX", required=False, default=100, par_type=int,
               docstr="Limit on the number of concurrent tasks in array")

        sf.par("NODESIZE", required=True, par_type=int,
               docstr="The number of cores per node defined by the system")

        sf.par("LSFARGS", required=False, default="", par_type=str,
               docstr="Any optional, additional LSG arguments that will be "
                      "passed to the LSF submit scripts")

    def submit(self, workflow):
        """
        Submits workflow
        """
        # Prepare 'bsub' arguments
        submit_call = " ".join([
            f"bsub",
            f"{PAR.LSFARGS}",
            f"-J {PAR.TITLE}",
            f"-o {self.output_log}.log",
            f"-e {self.error_log}.log",
            f"-n {PAR.NODESIZE}",
            f'-R "span[ptile={PAR.NODESIZE}"',
            f"-W {PAR.WALLTIME:d}:00",
            os.path.join(findpath("seisflows.system"), "wrappers", "submit"),
            PATH.OUTPUT
        ])

        super().submit(workflow, submit_call)

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
            f"bsub",
            f"{PAR.LSFARGS}",
            f"-J {PAR.TITLE}",
            f"-n {PAR.NPROC}",
            f'-R "span[ptile={PAR.NODESIZE}"',
            f"-W {PAR.TASKTIME:d}:00",
            f"-o {os.path.join(PATH.WORKDIR, 'output.logs', '%J_%I')}",
            f"[1-{PAR.NTASK}] % {PAR.NTASKMAX}",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}"
        ])

        stdout = subprocess.check_output(run_call, shell=True)

        # keep track of job ids
        jobs = self.job_id_list(stdout, PAR.NTASK)

        while True:
            # Wait seconds before checking status again
            time.sleep(30)
            self.timestamp()
            isdone, jobs = self.job_status(classname, method, jobs)
            if isdone:
                return

    def run_single(self, classname, method, *args, **kwargs):
        """ Runs task multiple times in embarrassingly parallel fasion

          Executes classname.method(*args, **kwargs) NTASK times, each time on
          NPROC cpu cores
        """
        # Checkpoint this individual method before proceeding
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        # Submit job array
        run_call = " ".join([
            f"bsub",
            f"{PAR.LSFARGS}",
            f"-J {PAR.TITLE}",
            f"-n {PAR.NPROC}",
            f'-R "span[ptile={PAR.NODESIZE}"',
            f"-W {PAR.TASKTIME:d}:00",
            f"-o {os.path.join(PATH.WORKDIR, 'output.logs', '%J')}",
            f"[1-1]",
            f"{os.path.join(findpath('seisflows.system'), 'wrappers', 'run')}",
            f"{PATH.OUTPUT}",
            f"{classname}",
            f"{method}",
            f"{PAR.ENVIRONS}"
        ])

        stdout = check_output(run_call, shell=True)

        # keep track of job ids
        jobs = self.job_id_list(stdout, ntask=1)

        while True:
            # Wait seconds before checking status again
            time.sleep(30)
            self.timestamp()
            isdone, jobs = self.job_status(classname, method, jobs)
            if isdone:
                return

    def job_id_list(self, stdout, ntask):
        """
        Parses job id list from sbatch standard output

        :type stdout: str
        :param stdout: the output of subprocess.check_output()
        :type ntask: int
        :param ntask: number of tasks currently running
        """
        job = stdout.split()[1].strip()[1:-1]
        if ntask == 1:
            return [job]
        else:
            number_jobs = range(1, PAR.NSRC + 1)
            return ["{job}[{}]".format(_) for _ in number_jobs]

    def job_status(self, classname, method, jobs):
        """
        Queries completion status of a single job

        :type job: str
        :param job: job id to query
        """
        job_finished = []
        for job in jobs:
            state = self._query(job)
            if state == "DONE":
                job_finished.append(True)
            else:
                job_finished.append(False)

            if state == "EXIT":
                print(msg.cli(f"LSF job {job} failed to execute "
                              f"{classname}.{method}.", header="error",
                              border="="))
                sys.exit(-1)

        isdone = all(job_finished)

        return isdone, jobs

    def _query(self, jobid):
        """
        Retrives job state from LSF database

        :type jobid: str
        :param jobid: job id to query LSF system about
        """
        # Write the job status output to a temporary file
        with open(os.path.join(PATH.SYSTEM, "job_status", "w")) as f:
            call('bjobs -a -d "{jobid}"', stdout=f)

        # Read the job status back from the text file
        with open(os.path.join(PATH.SYSTEM, "job_status", "r")) as f:
            lines = f.readlines()
            state = lines[1].split()[2].strip()

        return state

    def taskid(self):
        """
        Provides a unique identifier for each running task
        """
        return int(os.getenv('LSB_JOBINDEX')) - 1

    def timestamp(self):
        """
        Timestamp the current running job
        """
        with open(os.path.join(PATH.SYSTEM, "timestamps", "a")) as f:
            f.write(time.strftime("%H:%M:%S"))
            f.write("\n")

    def save_kwargs(self, classname, method, kwargs):
        """
        Save key word arguments as a pickle object.

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        """
        kwargspath = os.path.join(PATH.OUTPUT, "kwargs")
        kwargsfile = os.path.join(kwargspath, f"{classname}_{method}.p")

        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)
