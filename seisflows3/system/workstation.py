#!/usr/bin/env python3
"""
This is a subclass seisflows.system.workstation
Provides utilities for submitting jobs in serial on a single machine
"""
import os
import sys
import logging

from seisflows3.tools import unix, msg
from seisflows3.config import custom_import, SeisFlowsPathsParameters

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Workstation(custom_import("system", "base")):
    """
    Run tasks in a serial fashion on a single local machine
    """
    logger = logging.getLogger(__name__).getChild(__qualname__)

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
        self.setup()
        workflow.checkpoint()
        workflow.main()

    def run(self, classname, method, single=False, *args, **kwargs):
        """
        Executes task multiple times in serial.

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
        """
        unix.mkdir(PATH.SYSTEM)

        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        if single:
            os.environ["SEISFLOWS_TASKID"] = "0"
            func = getattr(__import__("seisflows_" + classname), method)
            func(**kwargs)
        else:
            for taskid in range(PAR.NTASK):
                os.environ["SEISFLOWS_TASKID"] = str(taskid)
                self.progress(taskid)
                func = getattr(__import__("seisflows_" + classname), method)
                func(**kwargs)

    def taskid(self):
        """
        Provides a unique identifier for each running task, which should be set
        by the 'run'' command.

        :rtype: int
        :return: returns the os environment variable SEISFLOWS_TASKID which is
            set by run() to label each of the currently
            running processes on the SYSTEM.
        """
        sftaskid = os.getenv("SEISFLOWS_TASKID")
        if sftaskid is None:
            print(msg.cli("system.taskid() environment variable not found. "
                          "Assuming DEBUG mode and returning taskid==0. "
                          "If not DEBUG mode, please check SYSTEM.run()",
                          header="warning", border="="))
            sftaskid = 0
        return int(sftaskid)

    def mpiexec(self):
        """
        Specifies MPI executable used to invoke solver

        .. note::
            For serial runs, MPIEXEC should be './' This is enforced in 
            tools.specfem.call_solver, which is the main function that uses
            PAR.MPIEXEC

        """
        return PAR.MPIEXEC

    def progress(self, taskid):
        """
        Provides status update by printing the current task being performed
        """
        if PAR.NTASK > 1:
            print(f"task {taskid + 1:02d} of {PAR.NTASK:02d}")
