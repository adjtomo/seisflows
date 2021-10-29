#!/usr/bin/env python
"""
This is a SeisFlows subclass which inherits attributes from a parent class
"""
import sys
from seisflows3.config import SeisFlowsPathsParameters, custom_import

# Required SeisFlows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

# The number of loaded modules depends on the module this Base class belongs to
system = sys.modules["seisflows_system"]
solver = sys.modules["seisflows_solver"]
optimize = sys.modules["seisflows_optimize"]
preprocess = sys.modules["seisflows_preprocess"]
postprocess = sys.modules["seisflows_postprocess"]


class Subclass(custom_import("MODULE NAME HERE", "PARENT CLASS NAME HERE")):
    """
    This is a template subclass
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.

        :rtype: seisflows.config.SeisFlowsPathsParameters
        :return: Paths and parameters that define the given class

        """
        # The super().required argument ensures that the sublcass inherits the
        # paths and parameters defined by its parent class
        sf = SeisFlowsPathsParameters(super().required)

        # > Additional or overloading paths and parameters can be set here

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths. The validate function ensures that all
        required paths and parameters are accounted for, and that all
        optional paths and parameters are set to user-defined or default values.
        """
        if validate:
            self.required.validate()
        # Validation only required by the lowest subclass, which will validate
        # all the paths and parameters from each of its parent classes
        super.check(validate=False)