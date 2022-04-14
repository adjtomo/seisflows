#!/usr/bin/env python3
"""
The System module provides the basic core utilities for interaction with compute
systems. The Base class must be overloaded by subclasses related to specific
compute system types (cluster vs. workstation) and even specific HPCs.
"""
import os
import sys
import pickle
import logging

from seisflows3.tools import unix, msg
from seisflows3.tools.wrappers import number_fid
from seisflows3.config import save, SeisFlowsPathsParameters, CFGPATHS



PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Base:
    """
    Abstract base class for the Systems module which controls interaction with
    compute systems such as HPC clusters.
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        self.output_log = os.path.join(PATH.WORKDIR, CFGPATHS.LOGFILE)
        self.error_log = os.path.join(PATH.WORKDIR, CFGPATHS.ERRLOGFILE)

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        sf.par("TITLE", required=False,
               default=os.path.basename(os.path.abspath(".")), par_type=str,
               docstr="The name used to submit jobs to the system, defaults "
                      "to the name of the working directory")

        sf.par("PRECHECK", required=False, par_type=list, default=["TITLE"],
               docstr="A list of parameters that will be displayed to stdout "
                      "before 'submit' or 'resume' is run. Useful for "
                      "manually reviewing important parameters prior to "
                      "system submission")

        sf.par("LOG_LEVEL", required=False, par_type=str, default="DEBUG",
               docstr="Verbosity output of SF3 logger. Available from least to "
                      "most verbosity: 'CRITICAL', 'WARNING', 'INFO', 'DEBUG'; "
                      "defaults to 'DEBUG'")

        sf.par("VERBOSE", required=False, default=True, par_type=bool,
               docstr="Level of verbosity provided to the output log. If True, "
                      "log statements will declare what module/class/function "
                      "they are being called from. Useful for debugging but "
                      "also very noisy.")

        # Define the Paths required by this module
        # note: PATH.WORKDIR has been set by the entry point seisflows.setup()
        sf.path("SCRATCH", required=False,
                default=os.path.join(PATH.WORKDIR, CFGPATHS.SCRATCHDIR),
                docstr="scratch path to hold temporary data during workflow")

        sf.path("OUTPUT", required=False,
                default=os.path.join(PATH.WORKDIR, CFGPATHS.OUTPUTDIR),
                docstr="directory to save workflow outputs to disk")

        sf.path("SYSTEM", required=False,
                default=os.path.join(PATH.WORKDIR, CFGPATHS.SCRATCHDIR,
                                     "system"),
                docstr="scratch path to hold any system related data")

        sf.path("LOCAL", required=False,
                docstr="path to local data to be used during workflow")

        sf.path("LOGFILE", required=False, default=self.output_log,
                docstr="the main output log file where all processes will "
                       "track their status")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        msg.check(type(self))

        if validate:
            self.required.validate()

        # !!! This system could be better, currently defining log file twice
        if self.output_log != PATH.LOGFILE:
            self.output_log = PATH.LOGFILE

    def setup(self):
        """
        Create the SeisFlows3 directory structure in preparation for a
        SeisFlows3 workflow. Ensure that if any config information is left over
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
        msg.setup(type(self))

        # Create scratch directories
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # Create output directories
        unix.mkdir(PATH.OUTPUT)
        log_files = os.path.join(PATH.WORKDIR, CFGPATHS.LOGDIR)
        unix.mkdir(log_files)

        # If resuming, move old log files to keep them out of the way. Number
        # in ascending order, so we don't end up overwriting things
        for src in [self.output_log, self.error_log, PATH.PAR_FILE]:
            i = 1
            if os.path.exists(src):
                dst = os.path.join(log_files, number_fid(src, i))
                while os.path.exists(dst):
                    i += 1
                    dst = os.path.join(log_files, number_fid(src, i))
                self.logger.debug(f"copying par/log file to: {dst}")
                unix.cp(src=src, dst=dst)

    def submit(self):
        """
        Main insertion point of SeisFlows3 onto the compute system.

        .. rubric::
            $ seisflows submit

        .. note::
            The expected behavior of the submit() function is to:
            1) run system setup, creating directory structure,
            2) execute workflow by submitting workflow.main()
        """
        raise NotImplementedError("Must be implemented by subclass")

    def run(self, classname, method, single=False, **kwargs):
        """
        Runs a task multiple times in am embarassingly parallel fashion

        .. note::
            The expected behavior of the run() function is to: submit N jobs to
            the system in parallel. For example, in a simulation step, run()
            submits N jobs to the compute system where N is the number of
            events requiring an adjoint simulation.

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
        :rtype: None
        :return: This function is not expected to return anything
        """
        raise NotImplementedError("Must be implemented by subclass")

    def taskid(self):
        """
        Provides a unique identifier for each running task. This is
        compute system specific.

        :rtype: int
        :return: this function is expected to return a unique numerical
            identifier.
        """
        raise NotImplementedError("Must be implemented by subclass")

    def checkpoint(self, path, classname, method, kwargs):
        """
        Writes the SeisFlows3 working environment to disk so that new tasks can
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
        self.logger.debug("checkpointing working environment to disk")

        argspath = os.path.join(path, "kwargs")
        argsfile = os.path.join(argspath, f"{classname}_{method}.p")
        unix.mkdir(argspath)
        with open(argsfile, "wb") as f:
            pickle.dump(kwargs, f)
        save()

