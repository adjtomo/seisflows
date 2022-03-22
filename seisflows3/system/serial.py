#!/usr/bin/env python
"""
This is a subclass seisflows.system.Serial
Provides utilities for submitting jobs in serial on a single machine, mostly
for testing purposes
"""
import os
import sys
import logging

from seisflows3.tools import unix, msg
from seisflows3.config import custom_import, SeisFlowsPathsParameters

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Serial(custom_import("system", "base")):
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

        sf.par("MPIEXEC", required=False, par_type=str,
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
        # Run setup to create necessary directory structure
        output_log, error_log = self.setup()
        workflow.checkpoint()

        # execute workflow
        workflow.main()

    def run(self, classname, method, hosts="all", *args, **kwargs):
        """
        Executes task multiple times in serial
        """
        unix.mkdir(PATH.SYSTEM)

        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        for taskid in range(PAR.NTASK):
            os.environ["SEISFLOWS_TASKID"] = str(taskid)
            self.progress(taskid)
            func = getattr(__import__("seisflows_" + classname), method)
            func(**kwargs)

    def run_single(self, classname, method, scale_tasktime=1, *args, **kwargs):
        """
        Runs task a single time, used for running serial tasks such as smoothing
        """
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)
        
        os.environ["SEISFLOWS_TASKID"] = "0"
        func = getattr(__import__("seisflows_" + classname), method)
        func(**kwargs)

    def taskid(self):
        """
        Provides a unique identifier for each running task, which should be set
        by the 'run' or 'run_single' command.

        :rtype: int
        :return: returns the os environment variable SEISFLOWS_TASKID which is
            set by run() or run_single() to label each of the currently
            running processes on the SYSTEM.
        """
        try:
            tid = int(os.environ["SEISFLOWS_TASKID"])
        except KeyError:
            # This should only return a KeyError if you're running in debug mode
            # and aren't assigned a TASKID by the OS. Return task id = 0 
            # i.e., mainsolver, so that user can efficiently run debug commands
            print(msg.cli("system.taskid() environment variable not found. "
                          "Assuming DEBUG mode and returning taskid==0. "
                          "If not DEBUG mode, please check "
                          "SYSTEM.run() or SYSTEM.run_single().",
                          header="warning", border="="))
            tid = 0

        return tid

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
