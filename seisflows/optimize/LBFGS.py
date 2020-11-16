#!/usr/bin/env python
"""
This is the custom class for an LBFGS optimization schema.
It supercedes the `seisflows.optimize.base` class
"""
import sys
import numpy as np

from seisflows.config import custom_import, SeisFlowsPathsParameters
from seisflows.plugins import optimize

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class LBFGS(custom_import("optimize", "base")):
    """
    The Limited memory BFGS algorithm
    Calls upon seisflows.plugin.optimize.LBFGS to accomplish LBFGS algorithm
    """
    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type LBFGS: Class
        :param LBFGS: plugin LBFGS class that controls the machinery of the
            L-BFGS optimization schema
        :type restarted: int
        :param restarted: a flag to let Seisflows know if the LBFGS algorithm
            has been restarted
        """
        super().__init__()
        self.LBFGS = None
        self.restarted = None

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        # Define the Parameters required by this module
        sf.par("LINESEARCH", required=False, default="Backtrack", par_type=str,
               docstr="Algorithm to use for line search, see "
                      "seisflows.plugins.line_search for available choices")

        sf.par("LBFGSMEM", required=False, default=3, par_type=int,
               docstr="Max number of previous gradients to retain "
                      "in local memory")

        sf.par("LBFGSMAX", required=False, par_type=int, default="inf",
               docstr="LBFGS periodic restart interval, between 1 and 'inf'")

        sf.par("LBFGSTHRESH", required=False, default=0., par_type=float,
               docstr="LBFGS angle restart threshold")

        return sf

    def check(self, validate=True):
        """
        Checks parameters, paths, and dependencies
        """
        super().check(validate=False)
        if validate:
            self.required.validate()

        assert(PAR.LINESEARCH.upper() == "BACKTRACK"), \
            "LBFGS requires a Backtracking line search"

    def setup(self):
        """
        Set up the LBFGS optimization schema
        """
        super().setup()
        self.LBFGS = getattr(optimize, "LBFGS")(path=PATH.OPTIMIZE,
                                                memory=PAR.LBFGSMEM,
                                                maxiter=PAR.LBFGSMAX,
                                                thresh=PAR.LBFGSTHRESH,
                                                precond=self.precond,
                                                verbose=PAR.VERBOSE)

    def compute_direction(self):
        """
        Overwrite the Base compute direction class
        """
        p_new, self.restarted = self.LBFGS()

        self.save("p_new", p_new)

    def restart(self):
        """
        Overwrite the Base restart class and include a restart of the LBFGS
        """
        super().restart()
        self.LBFGS.restart()

