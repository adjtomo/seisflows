#!/usr/bin/env python
"""
This is the custom class for an steepest descent optimization schema.
It supercedes the `seisflows.optimize.base` class
"""
import sys

from seisflows3.config import custom_import, SeisFlowsPathsParameters

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class SteepestDescent(custom_import("optimize", "base")):
    """
    Steepest descent method
    """
    def __init__(self):
        self.restarted = False

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        # Define the Parameters required by this module
        sf.overwrite("LINESEARCH", default="Bracket")

        return sf

    def check(self, validate=True):
        """
        Checks parameters, paths, and dependencies
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

    def setup(self):
        """
        Set up the steepest descent optimization schema
        """
        super().setup()

    def compute_direction(self):
        """
        Overwrite the Base compute direction class
        """
        super().compute_direction()

    def restart(self):
        """
        Steepest descent never requires restarts
        """
        pass

