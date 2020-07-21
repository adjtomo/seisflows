#!/usr/bin/env python
"""
This is a subclass seisflows.system.Serial
Provides utilities for submitting jobs in serial on a single machine, mostly
for testing purposes
"""
import os
import sys
import numpy as np

import os
from seisflows.tools import unix
from seisflows.config import custom_import
from seisflows.tools.err import ParameterError

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Serial(custom_import("system", "base")):
    """
    Run tasks in a serial fashion on a single local machine
    """
    def check(self):
        """
        Checks parameters and paths
        """
        # name of the job
        if "TITLE" not in PAR:
            setattr(PAR, "TITLE", os.path.basename(os.path.abspath(".")))

        # number of tasks
        if "NTASK" not in PAR:
            setattr(PAR, "NTASK", 1)

        # number of processers per task
        if "NPROC" not in PAR:
            setattr(PAR, "NPROC", 1)

        # how to invoke executables
        if "MPIEXEC" not in PAR:
            setattr(PAR, "MPIEXEC", "")

        # level of detail in output messages
        if "VERBOSE" not in PAR:
            setattr(PAR, "VERBOSE", 1)

        # where job was submitted
        if "WORKDIR" not in PATH:
            setattr(PATH, "WORKDIR", os.path.abspath("."))

        # where output files are written
        if "OUTPUT" not in PATH:
            setattr(PATH, "OUTPUT", os.path.join(PATH.WORKDIR, "output"))

        # where temporary files are written
        if "SCRATCH" not in PATH:
            setattr(PATH, "SCRATCH", os.path.join(PATH.WORKDIR, "scratch"))

        # where system files are written
        if "SYSTEM" not in PATH:
            setattr(PATH, "SYSTEM", os.path.join(PATH.SCRATCH, "system"))

        # optional local filesystem scratch path
        if "LOCAL" not in PATH:
            setattr(PATH, "LOCAL", None)

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
