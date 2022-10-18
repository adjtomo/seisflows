#!/usr/bin/env python3
"""
This is the subclass seisflows.system.lsf.Lsf

This class provides the core utilities interaction with HPC systems which run
using the Platform Load Sharing Facility (LSF) workload management platform.
"""
import os
import time
import subprocess
import sys
from seisflows import ROOT_DIR
from seisflows.system.cluster import Cluster


class Lsf(Cluster):
    """
    An interface through which to submit workflows, run tasks in serial or
    parallel, and perform other system functions.

    By hiding environment details behind a python interface layer, these
    classes provide a consistent command set across different computing
    environments.

    Intermediate files are written to a global scratch path self.path.SCRATCH,
    which must be accessible to all compute nodes.

    Optionally, users can provide a local scratch path self.path.LOCAL if each
    compute node has its own local filesystem.

    For important additional information, please see
    http://seisflows.readthedocs.org/en/latest/manual/
                                                manual.html#system-configuration
    """
    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        raise NotImplementedError("This module is still a work in progress")
        sys.exit(-1)

        super().__init__()

        self.logger.warning("system.LSF is underdeveloped and "
                            "will likely not work without significant testing "
                            "and source code edits")

        self.required.par(
            "MPIEXEC", required=False, default="mpiexec", par_type=str,
            docstr="Function used to invoke executables on the system. "
                   "For example 'srun' on SLURM systems, or './' on a "
                   "workstation. If left blank, will guess based on the "
                   "system."
        )
        # Define the Parameters required by this module
        self.required.par(
            "NTASKMAX", required=False, default=100, par_type=int,
            docstr="Limit on the number of concurrent tasks in array"
        )
        self.required.par(
            "NODESIZE", required=True, par_type=int,
            docstr="The number of cores per node defined by the system"
        )
        self.required.par(
            "LSFARGS", required=False, default="", par_type=str,
            docstr="Any optional, additional LSG arguments that will be "
                   "passed to the LSF submit scripts"
        )
        self.required.path(
            "LOCAL", required=False,
            docstr="path to local data to be used during workflow"
        )

    def submit(self, submit_call=None):
        """
        Submits workflow using 'bsub' arguments
        """
        if submit_call is None:
            submit_call = " ".join([
                f"bsub",
                f"{self.par.LSFARGS}",
                f"-J {self.par.TITLE}",
                f"-o {self.output_log}.log",
                f"-e {self.error_log}.log",
                f"-n {self.par.NODESIZE}",
                f'-R "span[ptile={self.par.NODESIZE}"',
                f"-W {self.par.WALLTIME:d}:00",
                f"{os.path.join(ROOT_DIR, 'scripts', 'submit')}",
                f"--output {self.path.OUTPUT}"
            ])

        super().submit(submit_call=submit_call)

    def run(self, classname, method, single=False, run_call=None, **kwargs):
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
        self.checkpoint(self.path.OUTPUT, classname, method, kwargs)

        # Submit job array
        run_call = " ".join([
            f"bsub",
            f"{self.par.LSFARGS}",
            f"-J {self.par.TITLE}",
            f"-n {self.par.NPROC}",
            f'-R "span[ptile={self.par.NODESIZE}"',
            f"-W {self.par.TASKTIME:d}:00",
            f"-o {os.path.join(self.path.WORKDIR, 'output.logs', '%J_%I')}",
            f"[1-{self.par.NTASK}] % {self.par.NTASKMAX}",
            f"{os.path.join(ROOT_DIR, 'scripts', 'run')}",
            f"--output {self.path.OUTPUT}"
            f"--classname {classname}",
            f"--funcname {method}",
            f"--environment {self.par.ENVIRONS or ''}"
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

        # keep track of job ids
        jobs = self._job_id_list(stdout, single)

        while True:
            # Wait seconds before checking status again
            time.sleep(5)
            isdone, jobs = self._job_status(classname, method, jobs)
            if isdone:
                return

    def taskid(self):
        """
        Provides a unique identifier for each running task
        """
        return int(os.getenv('LSB_JOBINDEX')) - 1

    def checkpoint(self, path, classname, method, kwargs):
        """Inherits from workflow.system.workstation.Workstation"""
        self.checkpoint(path=path, classname=classname, method=method,
                        kwargs=kwargs)

    def _check_job_status(self, job_ids):
        """
        Queries completion status of a single job

        TODO this function is mangled, needs to be rewritten

        :type job: str
        :param job: job id to query
        """
        job_finished = []
        for job_id in job_ids:
            state = self._query(job_id)
            if state == "DONE":
                job_finished.append(True)
            else:
                job_finished.append(False)

            if state == "EXIT":
                return job_id, "FAILED"

        isdone = all(job_finished)

        return None, "OKAY"

    def _job_id_list(self, stdout, single):
        """
        Parses job id list from LSF standard output

        :type stdout: str
        :param stdout: the output of subprocess.check_output()
        :type single: bool
        :param single: if running a single process job, returns a list of length
            1 with a single job id, else returns a list of length self.par.NTASK
            for all arrayed jobs
        :rtype: list
        :return: a list of array jobs that should be currently running
        """
        job = stdout.split()[1].strip()[1:-1]
        if single:
            return [job]
        else:
            number_jobs = range(1, self.par.NSRC + 1)
            return ["{job}[{}]".format(_) for _ in number_jobs]

    def _query(self, jobid):
        """
        Retrives job state from LSF database

        :type jobid: str
        :param jobid: job id to query LSF system about
        """
        # Write the job status output to a temporary file
        with open(os.path.join(self.path.SYSTEM, "job_status", "w")) as f:
            subprocess.call('bjobs -a -d "{jobid}"', stdout=f)

        # Read the job status back from the text file
        with open(os.path.join(self.path.SYSTEM, "job_status", "r")) as f:
            lines = f.readlines()
            state = lines[1].split()[2].strip()

        return state

    # def save_kwargs(self, classname, method, kwargs):
    #     """
    #     Save key word arguments as a pickle object.
    #
    #     :type classname: str
    #     :param classname: the class to run
    #     :type method: str
    #     :param method: the method from the given `classname` to run
    #     """
    #     kwargspath = os.path.join(self.path.OUTPUT, "kwargs")
    #     kwargsfile = os.path.join(kwargspath, f"{classname}_{method}.p")
    #
    #     unix.mkdir(kwargspath)
    #     saveobj(kwargsfile, kwargs)


def _modify_run_call_single_proc(run_call):
    """

    """
    raise NotImplementedError("This function needs to be written")
