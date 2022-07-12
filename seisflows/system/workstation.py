#!/usr/bin/env python3
"""
This is a subclass seisflows.system.workstation
Provides utilities for submitting jobs in serial on a single machine
"""
import os
import sys
import pickle
from contextlib import redirect_stdout

from seisflows import logger
from seisflows.config import CFGPATHS, save
from seisflows.tools import msg, unix
from seisflows.tools.wrappers import number_fid


class Workstation:
    """
    Run tasks in a serial fashion on a single local machine. Also serves as the
    Base System module, upon which all other System classes should be built.
    """
    def __init__(self, title=None, mpiexec=None, ntask=1, nproc=1,
                 log_level="DEBUG", verbose=False, path_output=None,
                 path_system=None, path_output_log=None, path_error_log=None,
                 path_log_files=None, path_par_file=None, **kwargs):
        """
        Instantiate the Workstation base class

        :type title: str
        :param title: The name used to submit jobs to the system, defaults
            to the name of the current working directory
        :type mpiexec: str
        :param mpiexec: Function used to invoke executables on the system.
            For example 'srun' on SLURM systems. If None this will default to
            './' for calling executables.
        :type ntask: int
        :param ntask: number of individual tasks/events to run during workflow
        :type nproc: int
        :param nproc: number of processors to use for each simulation
        :type log_level: str
        :param log_level: logger level to pass to logging module.
            Available: 'debug', 'info', 'warning'
        :type verbose: bool
        :param verbose: if True, formats the log messages to include the file
            name, line number and message type. Useful for debugging but
            also very verbose
        :type path_output: str
        :param path_output: path to save files permanently to disk
        :type path_system: str
        :param path_system: scratch path to save any system related files
        """
        self.title = title
        self.mpiexec = mpiexec
        self.ntask = ntask
        self.nproc = nproc
        self.log_level = log_level
        self.verbose = verbose

        # Define internal path system
        self.path = path_system or \
                    os.path.join(os.getcwd(), "scratch", "system")
        self.path_par_file = path_par_file or \
                             os.path.join(os.getcwd(), "parameters.yaml")
        self.path_output = path_output

        # Define where to write logs
        self.path_log_files = path_log_files or \
                              os.path.join(os.getcwd(), "logs")
        self.path_output_log = path_output_log or \
                               os.path.join(os.getcwd(), "sfoutput.log")
        self.path_error_log = path_error_log or \
                              os.path.join(os.getcwd(), "sferror.log")


    def check(self):
        """
        Checks parameters and paths
        """
        pass

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
        for path in [self.path, self.path_output, self.path_log_files]:
            unix.mkdir(path)

        # If resuming, move old log files to keep them out of the way. Number
        # in ascending order, so we don't end up overwriting things
        for src in [self.path_output_log, self.path_error_log,
                    self.path_par_file]:
            i = 1
            if os.path.exists(src):
                dst = os.path.join(self.path_log_files, number_fid(src, i))
                while os.path.exists(dst):
                    i += 1
                    dst = os.path.join(self.path_log_files, number_fid(src, i))
                logger.debug(f"copying par/log file to: {dst}")
                unix.cp(src=src, dst=dst)

    def submit(self, workflow, submit_call=None):
        """
        Submits the main workflow job as a serial job submitted directly to
        the compute node that is running the master job

        TO DO fix this

        :type submit_call: str or None
        :param submit_call: the command line workload manager call to be run by
            subprocess. This is only needed for overriding classes, it has no
            effect on the Workstation class
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
        self.save_kwargs_to_disk(self.path_output, classname, method, kwargs)

        # Allows dynamic retrieval of any function from within package, e.g.,
        # <bound method Base.eval_func of <seisflows.solver.specfem2d...
        class_module = sys.modules[f"seisflows_{classname}"]
        function = getattr(class_module, method)

        if single:
            ntasks = 1
        else:
            ntasks = self.ntask

        for taskid in range(ntasks):
            # os environment variables can only be strings, these need to be
            # converted back to integers by system.taskid()
            os.environ["SEISFLOWS_TASKID"] = str(taskid)

            # Make sure that we're creating new log files EACH time we run()
            idx = 0
            while True:
                log_file = os.path.join(self.path_log_files,
                                        f"{idx:0>4}_{taskid:0>2}.log")
                if os.path.exists(log_file):
                    idx += 1
                else:
                    break

            if taskid == 0:
                logger.info(f"running task {classname}.{method} "
                            f"{self.ntask} times")

            # Redirect output to a log file to mimic cluster runs where 'run'
            # task output logs are sent to different files
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

    def save_kwargs_to_disk(self, path, classname, method, kwargs):
        """
        Writes keyword arguments for a given method to disk

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

        save(path=self.path_output)