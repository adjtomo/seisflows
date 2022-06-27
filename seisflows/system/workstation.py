#!/usr/bin/env python3
"""
This is a subclass seisflows.system.workstation
Provides utilities for submitting jobs in serial on a single machine
"""
import os
import sys
import logging
from contextlib import redirect_stdout

from seisflows.tools import msg
from seisflows.core import SeisFlowsPathsParameters
from seisflows.config import custom_import, CFGPATHS

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

        sf.par("MPIEXEC", required=False, default=None, par_type=str,
               docstr="Function used to invoke executables on the system. "
                      "For example 'srun' on SLURM systems, or './' on a "
                      "workstation. If left blank, will guess based on the "
                      "system.")

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

    def submit(self):
        """
        Submits the main workflow job
        """
        self.setup()
        workflow = sys.modules["seisflows_workflow"]
        workflow.checkpoint()
        workflow.main()

    def run(self, classname, method, single=False, **kwargs):
        """
        Executes task multiple times in serial.

        .. note::
            kwargs will be passed to the underlying `method` that is called

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
        self.checkpoint(PATH.OUTPUT, classname, method, kwargs)

        # Allows dynamic retrieval of any function from within package, e.g.,
        # <bound method Base.eval_func of <seisflows.solver.specfem2d...
        class_module = sys.modules[f"seisflows_{classname}"]
        function = getattr(class_module, method)
        log_path = os.path.join(PATH.WORKDIR, CFGPATHS.LOGDIR)

        if single:
            ntasks = 1
        else:
            ntasks = PAR.NTASK

        for taskid in range(ntasks):
            # os environment variables can only be strings, these need to be
            # converted back to integers by system.taskid()
            os.environ["SEISFLOWS_TASKID"] = str(taskid)

            # Make sure that we're creating new log files EACH time we run()
            idx = 0
            while True:
                log_file = os.path.join(log_path, f"{idx:0>4}_{taskid:0>2}.log")
                if os.path.exists(log_file):
                    idx += 1
                else:
                    break

            if taskid == 0:
                self.logger.info(f"running task {classname}.{method} "
                                 f"{PAR.NTASK} times")
            # Redirect output to a log file to mimic cluster runs
            with open(log_file, "w") as f:
                with redirect_stdout(f):
                    function(**kwargs)

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
