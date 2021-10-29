#!/usr/bin/env python
"""
This is the base class seisflows.system.Base
This class provides the core utilities interaction with HPC systems which must
be overloaded by subclasses
"""
import os
import sys
from glob import glob
from seisflows3.tools import unix
from seisflows3.config import save, saveobj, SeisFlowsPathsParameters


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class Base:
    """
    Abstract base class for the Systems module which controls interaction with
    compute systems such as HPC clusters.
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        # Define the Parameters required by this module
        sf.par("TITLE", required=False,
               default=os.path.basename(os.path.abspath(".")), par_type=str,
               docstr="The name used to submit jobs to the system, defaults "
                      "to the name of the working directory")

        sf.par("WALLTIME", required=True, par_type=float,
               docstr="Maximum job time in minutes for main SeisFlows job")

        sf.par("TASKTIME", required=True, par_type=float,
               docstr="Maximum job time in minutes for each SeisFlows task")

        sf.par("NTASK", required=True, par_type=int,
               docstr="Number of separate, individual tasks. Also equal to "
                      "the number of desired sources in workflow")

        sf.par("NPROC", required=True, par_type=int,
               docstr="Number of processor to use for each simulation")

        sf.par("PRECHECK", required=False, par_type=list,
               default=["TITLE", "BEGIN", "END", "WALLTIME"],
               docstr="A list of parameters that will be displayed to stdout "
                      "before 'submit' or 'resume' is run. Useful for "
                      "manually reviewing important parameters prior to "
                      "system submission")

        # Define the Paths required by this module
        # note: PATH.WORKDIR has been set by the entry point seisflows.setup()
        sf.path("SCRATCH", required=False,
                default=os.path.join(PATH.WORKDIR, "scratch"),
                docstr="scratch path to hold temporary data during workflow")

        sf.path("OUTPUT", required=False,
                default=os.path.join(PATH.WORKDIR, "output"),
                docstr="directory to save workflow outputs to disk")

        sf.path("SYSTEM", required=False,
                default=os.path.join(PATH.WORKDIR, "scratch", "system"),
                docstr="scratch path to hold any system related data")


        sf.path("LOCAL", required=False,
                docstr="path to local data to be used during workflow")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        if validate:
            self.required.validate()

    def setup(self):
        """
        Create the SeisFlows directory structure in preparation for a
        SeisFlows workflow. Ensure that if any config information is left over
        from a previous workflow, that these files are not overwritten by
        the new workflow. Should be called by submit()
        """
        # Create scratch directories
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # Create output directories
        unix.mkdir(PATH.OUTPUT)
        unix.mkdir(os.path.join(PATH.WORKDIR, "output.logs"))
        unix.mkdir(os.path.join(PATH.WORKDIR, "logs.old"))

        # If a scratch directory is made outside the working directory, make
        # sure its accessible from within the working directory
        if not os.path.exists("./scratch"):
            unix.ln(PATH.SCRATCH, os.path.join(PATH.WORKDIR, "scratch"))

        output_log = os.path.join(PATH.WORKDIR, "output")
        error_log = os.path.join(PATH.WORKDIR, "error")

        # If resuming, move old log files to keep them out of the way
        for log in [output_log, error_log]:
            unix.mv(src=glob(os.path.join(f"{log}*.log")),
                    dst=os.path.join(PATH.WORKDIR, "logs.old")
                    )

        # Copy the parameter.yaml file into the log directoroy
        par_copy = f"parameters_{PAR.BEGIN}-{PAR.END}.yaml"
        unix.cp(src="parameters.yaml",
                dst=os.path.join(PATH.WORKDIR, "logs.old", par_copy)
                )

        return output_log, error_log

    def submit(self):
        """
        Submits workflow
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def run(self, classname, method, *args, **kwargs):
        """
        Runs task multiple times
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def run_single(self, classname, method, *args, **kwargs):
        """
        Runs task a single time
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def taskid(self):
        """
        Provides a unique identifier for each running task
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def checkpoint(self, path, classname, method, args, kwargs):
        """
        Writes information to disk so tasks can be executed remotely

        :type path: str
        :param path: path the kwargs
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
        saveobj(argsfile, kwargs)
        save()

