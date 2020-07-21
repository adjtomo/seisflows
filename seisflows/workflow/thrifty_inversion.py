#!/usr/bin/env python
"""
This is the subclass seisflows.workflow.thrifty_pyatoa.

This is a Seisflows subclass, it controls the main workflow.

This subclass inherits from the seisflows.workflow.inversion_pyatoa class
It allows the inversion to skip the costly intialization step if the final
forward simulations from the previous iteration can be used in the current one.

This is a direct copy of `ThriftyInversion`, but custom imports from
InversionPyatoa rather than Inversion
"""
import sys

from seisflows.tools import unix
from seisflows.config import custom_import

# Seisflows Configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

optimize = sys.modules['seisflows_optimize']


class ThriftyInversion(custom_import("workflow", "inversion")):
    """
    Thrifty inversion subclass for InversionPyatoa
    """
    def __init__(self):
        """
        :type status: int
        :param status: the current status of the inversion.
            0: First iteration, restart, or other conditions means inversion
               must default to normal behavior
            1: A well-scaled inversion can skip the function evaluation of the
               next iteration by using the previous iteration.
        """
        self.status = 0

    def initialize(self):
        """
        If line search can be carried over, skip initialization step
        Or if manually starting a new run, start with normal inversion init
        """
        if (self.status == 0) or (optimize.iter == PAR.BEGIN):
            super(ThriftyInversion, self).initialize()
        else:
            print("THRIFTY INITIALIZE")

    def clean(self):
        """
        Determine if forward simulation from line search can be carried over
        """
        self.update_status()

        if self.status == 1:
            print("THRIFTY CLEAN")
            unix.rm(PATH.GRAD)
            unix.mv(PATH.FUNC, PATH.GRAD)
            unix.mkdir(PATH.FUNC)
        else:
            super(ThriftyInversion, self).clean()

    def update_status(self):
        """
        Determine if line search forward simulation can be carried over
        """
        print("THRIFTY STATUS")
        # Only works for backtracking line search
        if PAR.LINESEARCH != "Backtrack":
            print("\t Line search not 'Backtrack', cannot run thrifty")
            self.status = 0
        # May not work on first iteration
        elif optimize.iter == PAR.BEGIN:
            print("\t First iteration of workflow, defaulting to Inversion")
            self.status = 0
        # May not work following restart
        elif optimize.restarted:
            print("\t Optimization has been restarted, defaulting to Inversion")
            self.status = 0
        # May not work after resuming saved workflow
        elif optimize.iter == PAR.END:
            print("\t End of workflow, defaulting to Inversion")
            self.status = 0
        # May not work if using local filesystems
        elif PATH.LOCAL:
            print("\t Local filesystem, cannot run thrifty")
            self.status = 0
        # Otherwise, continue with thrifty inversion
        else:
            print("\t Continuing with Thrifty Inversion")
            self.status = 1





