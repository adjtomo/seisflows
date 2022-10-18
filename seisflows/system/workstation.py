#!/usr/bin/env python3
"""
The `workstation` class is the foundational `System` module in SeisFlows,
it provides utilities for submitting jobs in SERIAL on a small-scale machine,
e.g., a workstation or a laptop. All other `System` classes build on this class.
"""
import os
import sys
import subprocess
from contextlib import redirect_stdout

from seisflows import logger
from seisflows.tools import unix
from seisflows.tools.config import Dict, import_seisflows
from seisflows.tools.config import number_fid, set_task_id


class Workstation:
    """
    Workstation System
    ------------------
    Defines foundational structure for System module. When used standalone,
    runs solver tasks either in serial (if `nproc`==1; i.e., without MPI) or in
    parallel (if `nproc`>1; i.e., with MPI). All other tasks are run in serial.

    Parameters
    ----------
    :type ntask: int
    :param ntask: number of individual tasks/events to run during workflow.
        Must be <= the number of source files in `path_specfem_data`
    :type nproc: int
    :param nproc: number of processors to use for each simulation. Choose 1 for
        serial simulations, and `nproc`>1 for parallel simulations.
    :type mpiexec: str
    :param mpiexec: MPI executable on system. Defaults to 'mpirun -n ${NPROC}'
    :type log_level: str
    :param log_level: logger level to pass to logging module.
        Available: 'debug', 'info', 'warning', 'critical'
    :type verbose: bool
    :param verbose: if True, formats the log messages to include the file
        name, line number and message type. Useful for debugging but
        also very verbose.

    Paths
    -----
    :type path_output_log: str
    :param path_output_log: path to a text file used to store the outputs of
        the package wide logger, which are also written to stdout
    :type path_par_file: str
    :param path_par_file: path to parameter file which is used to instantiate
        the package
    :type path_log_files: str
    :param path_log_files: path to a directory where individual log files are
        saved whenever a number of parallel tasks are run on the system.
    ***
    """
    def __init__(self, ntask=1, nproc=1, mpiexec=None, log_level="DEBUG",
                 verbose=False, workdir=os.getcwd(), path_output=None,
                 path_system=None, path_par_file=None, path_output_log=None,
                 path_log_files=None, **kwargs):
        """
        Workstation System Class Parameters

        .. note::
            Paths listed here are shared with `workflow.forward` and so are not
            included in the class docstring.

        :type workdir: str
        :param workdir: working directory in which to look for data and store
            results. Defaults to current working directory
        :type path_output: str
        :param path_output: path to directory used for permanent storage on disk.
            Results and exported scratch files are saved here.
        :type path_system: str
        :param path_system: scratch path to save any system related files
        """
        self.ntask = ntask
        self.nproc = nproc
        self.mpiexec = mpiexec
        self.log_level = log_level.upper()
        self.verbose = verbose

        # Define internal path system
        self.path = Dict(
            workdir=workdir or os.getcwd(),
            scratch=path_system or os.path.join(workdir, "scratch", "system"),
            par_file=path_par_file or os.path.join(workdir, "parameters.yaml"),
            output=path_output or os.path.join(workdir, "output"),
            log_files=path_log_files or os.path.join(workdir, "logs"),
            output_log=path_output_log or os.path.join(workdir, "sflog.txt"),
        )
        self._acceptable_log_levels = ["CRITICAL", "WARNING", "INFO", "DEBUG"]

    def check(self):
        """
        Checks parameters and paths
        """
        assert(os.path.exists(self.path.par_file)), \
            f"parameter file does not exist but should"

        assert(self.ntask > 0), f"number of events/tasks `ntask` cannot be neg'"
        assert(self.log_level in self._acceptable_log_levels), \
            f"`system.log_level` must be in {self._acceptable_log_levels}"

        if self.nproc > 1:
            assert(self.mpiexec is not None), (
                f"Multi-core workflows (`nproc`>1) require an MPI executable " 
                f"`mpiexec`"
            )
        if self.mpiexec is not None:
            # Make user that `mpiexec` exists on system
            try:
                stdout = subprocess.run(f"which {self.mpiexec}", shell=True,
                                        check=True, text=True,
                                        stdout=subprocess.PIPE).stdout
            except subprocess.CalledProcessError:
                logger.critical(
                    f"MPI executable {self.mpiexec} was not found on system "
                    f"with cmd: `which {self.mpiexec}. Please check that your "
                    f"MPI module is loaded and accessible from the command line"
                )
                sys.exit(-1)
            logger.debug(f"MPI executable is located at: {stdout.strip()}")

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
        for src in [self.path.output_log, self.path.par_file]:
            i = 1
            if os.path.exists(src):
                dst = os.path.join(self.path.log_files, number_fid(src, i))
                while os.path.exists(dst):
                    i += 1
                    dst = os.path.join(self.path.log_files, number_fid(src, i))
                logger.debug(f"copying par/log file to: {dst}")
                unix.cp(src=src, dst=dst)

    def submit(self, workdir=None, parameter_file="parameters.yaml"):
        """
        Submits the main workflow job as a serial job submitted directly to
        the system that is running the master job

        :type workdir: str
        :param workdir: path to the current working directory
        :type parameter_file: str
        :param parameter_file: parameter file name used to instantiate the
            SeisFlows package
        """
        workflow = import_seisflows(workdir=workdir or self.path.workdir,
                                    parameter_file=parameter_file)
        workflow.check()
        workflow.setup()
        workflow.run()

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

        for task_id in range(ntasks):
            # Set Task ID for currently running process
            set_task_id(task_id)
            log_file = self._get_log_file(task_id)

            # Redirect output to a log file to mimic cluster runs where 'run'
            # task output logs are sent to different files
            with open(log_file, "w") as f:
                with redirect_stdout(f):
                    for func in funcs:
                        func(**kwargs)

    def _get_log_file(self, task_id):
        """
        To mimic clusters which assign job numbers to spawned processes, our
        on-system runs will also assign job numbers simply be incrementing the
        number on the log files on system.
        """
        idx = 1
        while True:
            log_file = os.path.join(self.path.log_files,
                                    f"{idx:0>4}_{task_id:0>2}.log")
            if os.path.exists(log_file):
                idx += 1
            else:
                return log_file
