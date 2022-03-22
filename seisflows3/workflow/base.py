#!/usr/bin/env python
"""
This is the Base class for seisflows.workflow.
It contains mandatory functions that must be called by subclasses
"""
import logging
from seisflows3.config import save, SeisFlowsPathsParameters


class Base:
    """
    Workflow abstract base class
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        pass

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        sf.par("RESUME_FROM", required=False, par_type=str,
               docstr="Name of task to resume inversion from")

        sf.par("STOP_AFTER", required=False, par_type=str,
               docstr="Name of task to stop inversion after finishing")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths. Must be implemented by sub-class
        """
        pass

    def main(self):
        """
        Execution of a workflow is equal to stepping through workflow.main()
        """
        pass

    @staticmethod
    def checkpoint():
        """
        Writes information to disk so workflow can be resumed following a break
        """
        save()


