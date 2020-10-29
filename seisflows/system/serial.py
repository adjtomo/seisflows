#!/usr/bin/env python
"""
This is a subclass seisflows.system.Serial
Provides utilities for submitting jobs in serial on a single machine, mostly
for testing purposes
"""
import os
import sys

from seisflows.tools import unix
from seisflows.config import custom_import, SeisFlowsPathsParameters

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Serial(custom_import("system", "base")):
    """
    Run tasks in a serial fashion on a single local machine
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        # Define the Parameters required by this module
        sf.par("NTASK", required=False, default=1, par_type=int,
               docstr="Number of separate, individual tasks. Also equal to "
                      "the number of desired sources in workflow")

        sf.par("NPROC", required=False, default=1, par_type=int,
               docstr="Number of processor to use for each simulation")

        sf.par("MPIEXEC", required=False, default="", par_type=str,
               docstr="Function used to invoke parallel executables")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=False)
        if validate:
            self.required.validate()

    def submit(self, workflow):
        """
        Submits the main workflow job
        """
        # create scratch directories
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # create output directories
        unix.mkdir(PATH.OUTPUT)

        workflow.checkpoint()

        # execute workflow
        workflow.main()

    def run(self, classname, method, hosts="all", **kwargs):
        """
        Executes task multiple times in serial
        """
        unix.mkdir(PATH.SYSTEM)

        for taskid in range(PAR.NTASK):
            os.environ["SEISFLOWS_TASKID"] = str(taskid)
            self.progress(taskid)
            func = getattr(__import__("seisflows_" + classname), method)
            func(**kwargs)

    def run_single(self, classname, method, *args, **kwargs):
        """
        Runs task a single time
        """
        os.environ["SEISFLOWS_TASKID"] = "0"
        func = getattr(__import__("seisflows_" + classname), method)
        func(**kwargs)

    def taskid(self):
        """
        Provides a unique identifier for each running task
        """
        return int(os.environ["SEISFLOWS_TASKID"])

    def mpiexec(self):
        """
        Specifies MPI executable used to invoke solver
        """
        return PAR.MPIEXEC

    def progress(self, taskid):
        """
        Provides status update by printing the current task being performed
        """
        if PAR.NTASK > 1:
            print(f"task {taskid + 1:02d} of {PAR.NTASK:02d}")
