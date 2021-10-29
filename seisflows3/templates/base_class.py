#!/usr/bin/env python
"""
This is a SeisFlows Base class
"""
import os
import sys
from seisflows.config import SeisFlowsPathsParameters

# Required SeisFlows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

# The number of loaded modules depends on the module this Base class belongs to
system = sys.modules["seisflows_system"]
solver = sys.modules["seisflows_solver"]
optimize = sys.modules["seisflows_optimize"]
preprocess = sys.modules["seisflows_preprocess"]
postprocess = sys.modules["seisflows_postprocess"]


class Base:
    """
    This is a template Base class
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.

        :rtype: seisflows.config.SeisFlowsPathsParameters
        :return: Paths and parameters that define the given class

        """
        sf = SeisFlowsPathsParameters()

        # Define the Parameters required by this module
        sf.par("EXAMPLE_REQUIRED_PARAMETER", required=True, par_type=str,
               docstr="Required parameters do not need default values and will "
                      "need to be set by the user in the parameter file"
               )

        sf.par("EXAMPLE_OPTIONAL_PARAMETER", required=False, default=0,
               par_type=int, docstr="Optional parameters require a default "
                                    "value, if no default value is given, the "
                                    "parameter is set to None"
               )

        # Define the Paths required by this module
        sf.path("EXAMPLE_REQUIRED_PATH", required=True,
                docstr="Required paths to be set by user in parameter file")

        sf.path("EXAMPLE_OPTIONAL_PATH", required=False,
                default=os.path.join(PATH.SCRATCH, "example"),
                docstr="Optional paths require default values")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths. The validate function ensures that all
        required paths and parameters are accounted for, and that all
        optional paths and parameters are set to user-defined or default values.
        """
        if validate:
            self.required.validate()