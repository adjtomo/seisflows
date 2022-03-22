#!/usr/bin/env python
"""
This is a SeisFlows3 Test class which is used to test out the underlying
machinery before running an actual workflow. Contains simple functions used to
make sure that all parts of the package are working as expected.
"""
import os
import sys
import time
import logging
from seisflows3.tools import msg
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


class Test(custom_import("workflow", "base")):
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
        sf = SeisFlowsPathsParameters(super().required)

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

    def main(self, return_flow=False):
        """
        This controls the main testing workflow
        """
        FLOW = [self.test_system]
        if return_flow:
            return FLOW

        for func in FLOW:
            func()

    def test_function(self):
        """
        A simple function that can be called by system.run() or
        system.run_single()
        """
        print(f"Hello world, from taskid {system.taskid()}")

    def test_system(self):
        """
        This is an example test function which can take any number of args
        or kwargs. The base class is responsible for setting all of the
        necessary functions
        """
        system.run(classname="workflow", method="test_function")
        # Wait a bit for system to catch up
        time.sleep(3)
        system.run_single(classname="workflow", method="test_function")

    def 

