#!/usr/bin/env python
"""
This is the Base class for seisflows.workflow.
It contains mandatory functions that must be called by subclasses
"""
from seisflows.config import save, SeisFlowsPathsParameters


class Base:
    """
    Workflow abstract base class
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        return sf

    def check(self):
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


