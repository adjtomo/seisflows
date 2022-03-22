#!/usr/bin/env python
"""
This is the base class seisflows.system.cluster

This class provides the core utilities interaction with HPC systems which must
be overloaded by subclasses
"""
import os
import sys
import logging
from glob import glob

from seisflows3.tools import unix, msg
from seisflows3.tools.wrappers import call
from seisflows3.config import (custom_import, save, saveobj,
                               SeisFlowsPathsParameters, CFGPATHS)


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class Cluster(custom_import("system", "base")):
    """
    Abstract base class for the Systems module which controls interaction with
    compute systems such as HPC clusters.
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        # Define the Parameters required by this module
        sf.par("WALLTIME", required=True, par_type=float,
               docstr="Maximum job time in minutes for main SeisFlows3 job")

        sf.par("TASKTIME", required=True, par_type=float,
               docstr="Maximum job time in minutes for each SeisFlows3 task")

        sf.par("NTASK", required=True, par_type=int,
               docstr="Number of separate, individual tasks. Also equal to "
                      "the number of desired sources in workflow")

        sf.par("NPROC", required=True, par_type=int,
               docstr="Number of processor to use for each simulation")

        sf.par("ENVIRONS", required=False, default="", par_type=str,
               docstr="Optional environment variables to be provided in the"
                      "following format VAR1=var1,VAR2=var2...")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        msg.check(type(self))

        if validate:
            self.required.validate()

        super().check(validate=False)

    def submit(self, workflow, submit_call):
        """
        Main insertion point of SeisFlows3 onto the compute system.

        .. rubric::
            $ seisflows submit

        .. note::
            The expected behavior of the submit() function is to:
            1) run system setup, creating directory structure,
            2) execute workflow by submitting workflow.main()

        :type workflow: seisflows3.workflow
        :param workflow: an active seisflows3 workflow instance
        :type submit_call: str
        :param submit_call: the command line workload manager call to be run by
            subprocess. These need to be passed in by specific workload manager
            subclasses.
        """
        self.setup()
        workflow.checkpoint()
        call(submit_call)

    def run(self, classname, method, *args, **kwargs):
        """
        Runs a task multiple times in parallel

        .. note::
            The expected behavior of the run() function is to: submit N jobs to
            the system in parallel. For example, in a simulation step, run()
            submits N jobs to the compute system where N is the number of
            events requiring an adjoint simulation.

        :rtype: None
        :return: This function is not expected to return anything
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def run_single(self, classname, method, *args, **kwargs):
        """
        Runs a task a single time

        .. note::
            The expected behavior of the run_single() function is to submit ONE
            job to the compute system. This could be used for, submitting a job
            to smooth a model, which only needs to be done once.

        :rtype: None
        :return: This function is not expected to return anything
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def taskid(self):
        """
        Provides a unique identifier for each running task. This is
        compute system specific.

        :rtype: int
        :return: this function is expected to return a unique numerical
            identifier.
        """
        raise NotImplementedError('Must be implemented by subclass.')

    def checkpoint(self, path, classname, method, args, kwargs):
        """
        Writes the SeisFlows3 working environment to disk so that new tasks can
        be executed in a separate/new/restarted working environment.

        :type path: str
        :param path: path the kwargs
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
        saveobj(argsfile, kwargs)
        save()

