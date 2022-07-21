#!/usr/bin/env python3
"""
This is a subclass seisflows.system.workstation
Provides utilities for submitting jobs in serial on a single machine
"""
import os
from contextlib import redirect_stdout

from seisflows import logger
from seisflows.tools.core import Dict
from seisflows.tools import unix
from seisflows.tools.core import number_fid, get_task_id, set_task_id


class Workstation:
    """
    [system.workstation] runs tasks in serial on a local machine.

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
    def __init__(self, ntask=1, nproc=1, log_level="DEBUG", verbose=False,
                 workdir=os.getcwd(), path_output=None, path_system=None,
                 path_output_log=None, path_error_log=None, path_log_files=None,
                 path_par_file=None, **kwargs):
        """Workstation System Class Parameters"""
        self.ntask = ntask
        self.nproc = nproc
        self.log_level = log_level
        self.verbose = verbose

        # Define internal path system
        self.path = Dict(
            scratch=path_system or os.path.join(workdir, "scratch", "system"),
            par_file=path_par_file or os.path.join(workdir, "parameters.yaml"),
            output=path_output or os.path.join(workdir, "output"),
            log_files=path_log_files or os.path.join(workdir, "logs"),
            output_log=path_output_log or os.path.join(workdir, "sfoutput.log"),
            error_log=path_error_log or os.path.join(workdir, "sferror.log"),
        )

    def check(self):
        """
        Checks parameters and paths
        """
        assert(os.path.exists(self.path.par_file)), \
            f"parameter file does not exist but should"

    def taskid(self):
        """
        Provides a unique identifier for each running task, which should be set
        by the 'run'' command.

        :rtype: int
        :return: returns the os environment variable SEISFLOWS_TASKID which is
            set by run() to label each of the currently
            running processes on the SYSTEM.
        """
        return get_task_id()

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
        for path in [self.path.scratch, self.path.output, self.path.log_files]:
            unix.mkdir(path)

        # If resuming, move old log files to keep them out of the way. Number
        # in ascending order, so we don't end up overwriting things
        for src in [self.path.output_log, self.path.error_log,
                    self.path.par_file]:
            i = 1
            if os.path.exists(src):
                dst = os.path.join(self.path.log_files, number_fid(src, i))
                while os.path.exists(dst):
                    i += 1
                    dst = os.path.join(self.path.log_files, number_fid(src, i))
                logger.debug(f"copying par/log file to: {dst}")
                unix.cp(src=src, dst=dst)

    # def submit(self, workflow, submit_call=None):
    #     """
    #     Submits the main workflow job as a serial job submitted directly to
    #     the compute node that is running the master job
    #
    #     TO DO fix this
    #
    #     :type submit_call: str or None
    #     :param submit_call: the command line workload manager call to be run by
    #         subprocess. This is only needed for overriding classes, it has no
    #         effect on the Workstation class
    #     """
    #     self.setup()
    #     workflow = sys.modules["seisflows_workflow"]
    #     workflow.checkpoint()
    #     workflow.main()

    def run(self, funcs, single=False, **kwargs):
        """
        Executes task multiple times in serial.

        .. note::
            kwargs will be passed to the underlying `method` that is called

        :type funcs: list of methods
        :param funcs: a list of functions that should be run in order. All
            kwargs passed to run() will be passed into the functions.
        :type single: bool
        :param single: run a single-process, non-parallel task, such as
            smoothing the gradient, which only needs to be run by once.
            This will change how the job array and the number of tasks is
            defined, such that the job is submitted as a single-core job to
            the system.
        """
        if single:
            ntasks = 1
        else:
            ntasks = self.ntask

        for taskid in range(ntasks):
            # Set Task ID for currently running process
            set_task_id(taskid)

            # Make sure that we're creating new log files EACH time we run()
            idx = 0
            while True:
                log_file = os.path.join(self.path.log_files,
                                        f"{idx:0>4}_{taskid:0>2}.log")
                if os.path.exists(log_file):
                    idx += 1
                else:
                    break

            # Redirect output to a log file to mimic cluster runs where 'run'
            # task output logs are sent to different files
            with open(log_file, "w") as f:
                with redirect_stdout(f):
                    for func in funcs:
                        func(**kwargs)
