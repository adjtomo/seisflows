#!/usr/bin/env python3
"""
This is a subclass seisflows.system.workstation
Provides utilities for submitting jobs in serial on a single machine
"""
import os
import sys
import pickle
from contextlib import redirect_stdout

from seisflows.core import Base
from seisflows.config import CFGPATHS, save
from seisflows.tools import msg, unix
from seisflows.tools.wrappers import number_fid


class Workstation(Base):
    """
    Run tasks in a serial fashion on a single local machine. Also serves as the
    Base System module, upon which all other System classes should be built.
    """
    def __init__(self):
        """
        Instantiate the Workstation base class
        """
        super().__init__()

        self.output_log = os.path.join(self.path.WORKDIR, CFGPATHS.LOGFILE)
        self.error_log = os.path.join(self.path.WORKDIR, CFGPATHS.ERRLOGFILE)

        self.required.par(
            "TITLE", required=False,
            default=os.path.basename(os.path.abspath(".")), par_type=str,
            docstr="The name used to submit jobs to the system, defaults "
                   "to the name of the working directory"
        )
        self.required.par(
            "MPIEXEC", required=False, default=None, par_type=str,
            docstr="Function used to invoke executables on the system. "
                   "For example 'srun' on SLURM systems, or './' on a "
                   "workstation. If left blank, will guess based on the "
                   "system."
        )
        self.required.par(
            "NTASK", required=False, default=1, par_type=int,
            docstr="Number of separate, individual tasks. Also equal to "
                   "the number of desired sources in workflow"
        )
        self.required.par(
            "NPROC", required=False, default=1, par_type=int,
            docstr="Number of processor to use for each simulation"
        )
        self.required.par(
            "PRECHECK", required=False, par_type=list, default=["TITLE"],
            docstr="A list of parameters that will be displayed to stdout "
                   "before 'submit' or 'resume' is run. Useful for "
                   "manually reviewing important parameters prior to "
                   "system submission"
        )
        self.required.par(
            "LOG_LEVEL", required=False, par_type=str, default="DEBUG",
            docstr="Verbosity output of SF logger. Available from least to "
                   "most verbosity: 'CRITICAL', 'WARNING', 'INFO', 'DEBUG'; "
                   "defaults to 'DEBUG'"
        )
        self.required.par(
            "VERBOSE", required=False, default=False, par_type=bool,
            docstr="Level of verbosity provided to the output log. If True, "
                   "log statements will declare what module/class/function "
                   "they are being called from. Useful for debugging but "
                   "also very noisy."
        )
        # note: self.path.WORKDIR has been set by the entry point seisflows.setup()
        self.required.path(
            "SCRATCH", required=False,
            default=os.path.join(self.path.WORKDIR, "scratch"),
            docstr="scratch path to hold temporary data during workflow"
        )
        self.required.path(
            "OUTPUT", required=False,
            default=os.path.join(self.path.WORKDIR, "output"),
            docstr="directory to save workflow outputs to disk"
        )
        self.required.path(
            "SYSTEM", required=False,
            default=os.path.join(self.required.SCRATCH, "system"),
            docstr="scratch path to hold any system related data"
        )
        self.required.path(
            "LOGFILE", required=False, default=self.output_log,
            docstr="the main output log file where all processes will track "
                   "their status"
        )

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

        if self.output_log != self.path.LOGFILE:
            self.output_log = self.path.LOGFILE

    def setup(self):
        """
        Create the SeisFlows directory structure in preparation for a
        SeisFlows workflow. Ensure that if any config information is left over
        from a previous workflow, that these files are not overwritten by
        the new workflow. Should be called by submit()

        .. note::
            This function is expected to create dirs: SCRATCH, SYSTEM, OUTPUT
            and the following log files: output, error

        .. note::
            Logger is configured here as all workflows, independent of system,
            will be calling setup()

        :rtype: tuple of str
        :return: (path to output log, path to error log)
        """
        # Create scratch directories
        unix.mkdir(self.path.SCRATCH)
        unix.mkdir(self.path.SYSTEM)

        # Create output directories
        unix.mkdir(self.path.OUTPUT)
        log_files = os.path.join(self.path.WORKDIR, CFGPATHS.LOGDIR)
        unix.mkdir(log_files)

        # If resuming, move old log files to keep them out of the way. Number
        # in ascending order, so we don't end up overwriting things
        for src in [self.output_log, self.error_log, self.path.PAR_FILE]:
            i = 1
            if os.path.exists(src):
                dst = os.path.join(log_files, number_fid(src, i))
                while os.path.exists(dst):
                    i += 1
                    dst = os.path.join(log_files, number_fid(src, i))
                self.logger.debug(f"copying par/log file to: {dst}")
                unix.cp(src=src, dst=dst)

    def finalize(self):
        """Inherits from seisflows.core.Base"""
        super().finalize()

    def submit(self, submit_call=None):
        """
        Submits the main workflow job as a serial job submitted directly to
        the compute node that is running the master job

        :type submit_call: str or None
        :param submit_call: the command line workload manager call to be run by
            subprocess. This is only needed for overriding classes, it has no
            effect on the Workstation class
        """
        self.setup()
        workflow = self.module("workflow")
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
        self.checkpoint(self.path.OUTPUT, classname, method, kwargs)

        # Allows dynamic retrieval of any function from within package, e.g.,
        # <bound method Base.eval_func of <seisflows.solver.specfem2d...
        class_module = self.module(classname)
        function = getattr(class_module, method)
        log_path = os.path.join(self.path.WORKDIR, CFGPATHS.LOGDIR)

        if single:
            ntasks = 1
        else:
            ntasks = self.par.NTASK

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
                                 f"{self.par.NTASK} times")
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

    def checkpoint(self, path, classname, method, kwargs):
        """
        Writes the SeisFlows working environment to disk so that new tasks can
        be executed in a separate/new/restarted working environment.

        :type path: str
        :param path: path to save the checkpointed pickle files to
        :type classname: str
        :param classname: name of the class to save
        :type method: str
        :param method: the specific function to be checkpointed
        :type kwargs: dict
        :param kwargs: dictionary to pass to object saving
        """
        argspath = os.path.join(path, "kwargs")
        argsfile = os.path.join(argspath, f"{classname}_{method}.p")

        unix.mkdir(argspath)
        with open(argsfile, "wb") as f:
            pickle.dump(kwargs, f)
        save(path=self.path.OUTPUT)