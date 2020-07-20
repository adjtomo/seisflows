#!/usr/bin/env python
"""
This is the Base class for seisflows.workflow.
It contains mandatory functions that must be called by subclasses
"""
import time
from seisflows.config import save


class Base(object):
    """
    Workflow abstract base class
    """
    class Decorators:
        """
        A class of Decorator functions that can be general to various Workflows
        """
        @classmethod
        def stopwatch(cls, function):
            """
            Print out the time it took to run a given function. Useful for log
            statements 

            :type method: str
            :param method: interact with the stopwatch:
                "set" to start the timer, "time" to return time since `set` time
                or None to return the current time
            :rtype: time.time
            :return: a time object
            """
            start = time.time()
            function()
            print(f"{(time.time() - start) / 60.:.2f}m elapsed")

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


