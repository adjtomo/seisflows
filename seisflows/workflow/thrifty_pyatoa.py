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


class ThriftyPyatoa(custom_import('workflow', 'inversion_pyatoa')):
    """
    Thrifty inversion subclass for InversionPyatoa
    """
    # Instanstiate the status attribute
    status = 0

    def check(self):
        """
        Checks parameters and paths
        """
        # Run Base class checks
        super(ThriftyPyatoa, self).check()

        # Signifiy if data-synth. or synth.-synth. case
        if "FORCE_THRIFTY" not in PAR:
            setattr(PAR, "FORCE_THRIFTY", False)

    def initialize(self):
        """
        If line search can be carried over, skip initialization step
        Or if manually starting a new run, start with normal inversion init
        """
        if (self.status == 0) or (optimize.iter == PAR.BEGIN):
            super(ThriftyPyatoa, self).initialize()
        else:
            print("THRIFTY INITIALIZE")

    def clean(self, status=None):
        """
        Determine if forward simulation from line search can be carried over
        Allow manual control over status, so that the User can force a clean    
        Useful, e.g. if you force a thrifty inversion at PAR.END but end up 
        changing some parameters.
    
        :type status: int
        :param status: 0 or 1, 0 means default Inversion, clean scratch 
                               1 means Thrifty Inversion, do not clean 
        """
        if status is None:
            self.update_status()
        else:
            self.status = status
        
        if self.status == 1:
            print("THRIFTY CLEAN")
            unix.rm(PATH.GRAD)
            unix.mv(PATH.FUNC, PATH.GRAD)
            unix.mkdir(PATH.FUNC)
        else:
            super(ThriftyPyatoa, self).clean()

    def update_status(self):
        """
        Determine if line search forward simulation can be carried over
        """
        print("THRIFTY STATUS")
        # Only works for backtracking line search
        if PAR.LINESEARCH != "Backtrack":
            print("\t Line search not 'Backtrack', cannot run ThriftyInversion")
            self.status = 0
        # May not work on first iteration, allow force if no parameters changed
        elif optimize.iter == 1:
            print("\t L-BFGS requires 2 gradient evaluations for scaling, "
                  "defaulting to Inversion")
            self.status = 0
        # May not work following restart
        elif optimize.restarted:
            print("\t Optimization has been restarted, defaulting to Inversion")
            self.status = 0
        # May not work after resuming saved workflow
        elif optimize.iter == PAR.END and not PAR.FORCE_THRIFTY:
            print("\t End of workflow, defaulting to Inversion")
            self.status = 0
        # May not work if using local filesystems
        elif PATH.LOCAL:
            print("\t Local filesystem, cannot run ThriftyInversion")
            self.status = 0
        # Otherwise, continue with thrifty inversion
        else:
            print("\t Continuing with ThriftyInversion")
            self.status = 1



