#!/usr/bin/env python
"""
This is the Base class for seisflows.workflow.
It contains mandatory functions that must be called by subclasses
"""
from seisflows.config import save


class Base(object):
    """
    Workflow abstract base class
    """

    def check(self):
        """
        Checks parameters and paths
        """
        raise NotImplementedError("Must be implemented by subclass.")

    def main(self):
        """
        Main routine

        Execution of a workflow is equivalent to stepping through workflow.main
        """
        raise NotImplementedError("Must be implemented by subclass.")

    def checkpoint(self):
        """
        Writes information to disk so workflow can be resumed following a break
        """
        save()


