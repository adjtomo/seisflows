#!/usr/bin/env python
"""
This is a SeisFlows Base class
"""
import os
import sys
import logging
from seisflows3.tools import msg
from seisflows3.config import SeisFlowsPathsParameters

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
    # Class-specific logger accessed using self.logger
    # When this logger is called, e.g., self.logger.info("text"), the logging
    # package will know exactly which module, class and function the log
    # statement has been sent from, extraordinarily helpful for debugging.
    logger = logging.getLogger(__name__).getChild(__qualname__)

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

        :type validate: bool
        :param validate: set required paths and parameters into sys.modules
        """
        # Call to output log statement identifying this specific module + class
        msg.check(type(self))

        # The validate statement is used internally to set required paths
        # and parameters into sys.modules. Default values are stored for
        # optional terms
        if validate:
            self.required.validate()

    def test(self, *args, **kwargs):
        """
        This is an example test function which can take any number of args
        or kwargs. The base class is responsible for setting all of the
        necessary functions
        """
        #
        super.test()
        # Multiple logging levels determine how verbose the module will be
        self.logger.info("important log statement goes here")
        self.logger.debug("debugging log statement goes here")
        self.logger.warning("warnings can be passed here")
