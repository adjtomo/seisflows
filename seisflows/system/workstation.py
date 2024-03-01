#!/usr/bin/env python3
"""
The `workstation` class is the foundational `System` module in SeisFlows,
it provides utilities for submitting jobs in SERIAL on a small-scale machine,
e.g., a workstation or a laptop. All other `System` classes build on this class.
"""
import os
import sys
import subprocess
import time
import numpy as np
from contextlib import redirect_stdout

from seisflows import logger
from seisflows.tools import unix, msg
from seisflows.tools.config import Dict, import_seisflows
from seisflows.tools.config import copy_file, set_task_id


class Workstation:
    """
    Workstation System [System Base]
    --------------------------------
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
    :type tasktime: float
    :param tasktime: maximum job time in units minutes for each job spawned by
        the SeisFlows master job during a workflow. These include, e.g.,
        running the forward solver, adjoint solver, smoother, kernel combiner.
        All spawned tasks receive the same task time. Fractions of minutes
        acceptable. If set as `None`, no tasktime will be enforced.
    :type mpiexec: str
    :param mpiexec: MPI executable on system. Defaults to 'mpirun -n ${NPROC}'
    :type array: str
    :param array: for `ntask` > 1, determine which tasks to submit to run. By
        default (NoneType) this submits all task IDs [0:ntask), or for single
        runs, submits only the first task ID, 0. However, for debugging or
        manual control purposes, Users may input a string of task IDs that they
        would like to run. Follows formatting of SLURM array directive
        (https://slurm.schedmd.com/job_array.html), which is, for example:
        1,2,3-8:2,10 -> 1,2,3,5,7,10
        where '-' denotes a range (inclusive), and ':' denotes an optional step.
        If ':' step is not given for a range, then step defaults to 1.
    :type rerun: int
    :param rerun: [EXPERIMENTAL FEATURE] attempt to re-run failed tasks or 
        array tasks submitted with `run`. Collects information about failed 
        jobs (or array jobs) after a failure, and re-submits with `run`. 
        `rerun` is an integer defining how many times the User wants System to
        try and rerun before failing the entire job. If 0 (default), a single
        task failure will cause main job failure.
    :type log_level: str
    :param log_level: logger level to pass to logging module.
        Available: 'debug', 'info', 'warning', 'critical'
    :type verbose: bool
    :param verbose: if True, formats the log messages to include the file
        name and line number of the log message in the source code, as well as 
        the message and message type. Useful for debugging but also very verbose
        so not recommended for production runs.

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
    def __init__(self, ntask=1, nproc=1, tasktime=1, mpiexec=None, 
                 array=None, rerun=0, log_level="DEBUG", verbose=False, 
                 workdir=os.getcwd(), path_output=None, path_system=None, 
                 path_par_file=None, path_output_log=None, path_log_files=None, 
                 **kwargs):
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
        self.tasktime = tasktime
        self.rerun = rerun
        self.mpiexec = mpiexec
        self.array = array
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
                subprocess.run(f"which {self.mpiexec}", shell=True, check=True,
                               text=True, stdout=subprocess.PIPE).stdout
            except subprocess.CalledProcessError:
                logger.critical(
                    f"MPI executable '{self.mpiexec}' was not found on system "
                    f"with cmd: `which {self.mpiexec}`. Please check that your "
                    f"MPI module is loaded and accessible from the command line"
                )
                sys.exit(-1)
        # Check formatting of array parameters if provided
        if self.array is not None:
            try:
                self.task_ids()
            except Exception as e:
                logger.critical(f"`array` argument can not be parsed by System "
                                f"module. Please check error message: {e}")
                sys.exit(-1)
            logger.info(f"`system.array` == {self.array}")
    
        assert(isinstance(self.rerun, int)), f"`rerun` must be an int [0,inf)"
        assert(self.rerun >= 0), f"`rerun` must be in bounds [0, inf)"

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

    def finalize(self):
        """Tear down tasks for the end of an Inversion-based iteration"""
        unix.rm(self.path.scratch)
        unix.mkdir(self.path.scratch)

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
        # Copy log files if present to avoid overwriting
        for src in [self.path.output_log, self.path.par_file]:
            if os.path.exists(src) and os.path.exists(self.path.log_files):
                copy_file(src, copy_to=self.path.log_files)

        workflow = import_seisflows(workdir=workdir or self.path.workdir,
                                    parameter_file=parameter_file)
        workflow.check()
        workflow.setup()
        workflow.run()

    def run(self, funcs, single=False, tasktime=None, **kwargs):
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
        :type tasktime: float
        :param tasktime: Custom tasktime in units minutes for running the given 
            functions `funcs`. If not given, defaults to the System variable
            `tasktime`. If System `tasktime` is also None, defaults to no 
            tasktime (inifinty time). If tasks exceed the given `tasktime`, 
            the program will exit
        """
        if tasktime is None:
            tasktime = self.tasktime or np.inf  

        for task_id in self.task_ids(single):
            _start = time.time()
            # Set Task ID for currently running process
            set_task_id(task_id)
            log_file = self._get_log_file(task_id)

            # Redirect output to a log file to mimic cluster runs where 'run'
            # task output logs are sent to different files
            with open(log_file, "w") as f:
                with redirect_stdout(f):
                    for func in funcs:
                        func(**kwargs)
                        # Check timeout condition. This is imprecise as it 
                        # allows the function to finish first. But we assume 
                        # workstation approaches don't need precise timekeeping
                        if time.time() - _start > (tasktime * 60):  # seconds
                            logger.critical(msg.cli(
                                f"System Task {task_id} has exceeded the "
                                f"defined tasktime {tasktime}m. Please check "
                                f"log files and adjust `tasktime` if required", 
                                header="timeout error",  border="=")
                                )
                            sys.exit(-1)

    def task_ids(self, single=False):
        """
        Return a list of Task IDs (linked to each indiviudal source) to supply
        to the 'run' function. By default this returns a range of available
        tasks [0:ntask). See class docstring of parameter `array` for how to
        manually set task_ids to use for run call.

        :type single: bool
        :param single: If we only want to run a single process, this is will 
            default to TaskID == 0
        :rtype: list
        :return: a list of task IDs to be used by the `run` function
        """
        if single:
            task_ids = [0]
        else:
            if self.array is not None:
                task_ids = []
                parts = self.array.split(",")  # e.g., 1,3,5-9,10
                for part in parts:
                    # e.g., 5-9 -> 5,6,7,8,9
                    if "-" in part:
                        # e.g., 5-9:2 -> 5,7,9
                        if ":" in part:
                            part, step = part.split(":")
                        else:
                            step = 1
                        start, stop = part.split("-")
                        # Range is inclusive of both ends, default step is 1
                        task_ids += list(range(int(start), int(stop) + 1, 
                                               int(step)))
                    else:
                        task_ids.append(int(part))

            else:
                task_ids = list(range(0, self.ntask, 1))

        return task_ids

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
