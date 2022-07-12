#!/usr/bin/env python3
"""
The simplest simulation workflow you can run is a large number of forward
simulations to generate synthetics from a velocity model. Therefore the
Forward class represents the BASE workflow. All other workflows will build off
of the scaffolding defined by the Forward class.
"""
from seisflows.config import import_seisflows


class Base(object):
    """
    Defines the core Base object for all SeisFlows modules. All modules MUST
    inherit from the Base object to work properly. This Base class essentially
    dictates the required structure of a SeisFlows class.
    """
    def __init__(self):
        """
        SeisFlows instantiates its required parameters through the
        SeisFlowsPathsParameters class, which scaffolds a rigid framework of
        how parameters and paths should be defined by the program. This is
        then used to build the parameter file dynamically.
        """
        self.parameters = self.modules = import_seisflows()
        (self.system, self.preprocess, self.solver,
         self.postprocess, self.optimize) = self.modules

    def setup(self):
        """

        """

    @property
    def par(self):
        """
        Quick access SeisFlows parameters from sys.modules. Throws a warning
        if parameters have not been instantiated

        :rtype: Dict or None
        :return: Returns a Dictionary with instantiated parameters, or None if
            parameters have not been instantiated
        """
        return self.module("parameters")

    @property
    def path(self):
        """
        Quick access SeisFlows paths from sys.modules. Throws a warning
        if paths have not been instantiated

        :rtype: Dict or None
        :return: Returns a Dictionary with instantiated paths, or None if
            parameters have not been instantiated
        """
        return self.module("paths")

    def check(self, validate=True):
        """
        General check() function for each module to check the validity of the
        user-input parameters and paths
        """
        if validate:
            self.required.validate()

        # Example of a check statement
        # assert(self.par.PARAMETER == example_value), f"Parameter != example"

    def setup(self):
        """
        A placeholder function for any initialization or setup tasks that
        need to be run once at the beginning of any workflow.
        """
        return

    def finalize(self):
        """
        A placeholder function for any finalization or tear-down tasks that
        need to be run at the end of any iteration or workflow.
        """
        return