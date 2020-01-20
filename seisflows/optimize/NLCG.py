#!/usr/bin/env python
"""
This is the custom class for an NLCG optimization schema.
It supercedes the `seisflows.optimize.base` class
"""
import sys
import numpy as np

from seisflows.config import custom_import
from seisflows.plugins import optimize

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class NLCG(custom_import("optimize", "Base")):
    """
    Nonlinear conjugate gradient method
    """

    def check(self):
        """
        Checks parameters, paths, and dependencies
        """
        # Line search algorithm
        if "LINESEARCH" not in PAR:
            setattr(PAR, "LINESEARCH", "Bracket")

        # NLCG periodic restart interval
        if "NLCGMAX" not in PAR:
            setattr(PAR, "NLCGMAX", np.inf)

        # NLCG conjugacy restart threshold
        if "NLCGTHRESH" not in PAR:
            setattr(PAR, "NLCGTHRESH", np.inf)

        super(NLCG, self).check()

    def setup(self):
        """
        Set up the NLCG optimization schema
        """
        super(NLCG, self).setup()
        self.NLCG = getattr(optimize, 'NLCG')(path=PATH.OPTIMIZE,
                                              maxiter=PAR.NLCGMAX,
                                              thresh=PAR.NLCGTHRESH,
                                              precond=self.precond)

    def compute_direction(self):
        """
        Overwrite the Base compute direction class
        """
        # g_new = self.load('g_new')
        p_new, self.restarted = self.NLCG()
        self.save('p_new', p_new)

    def restart(self):
        """
        Overwrite the Base restart class and include a restart of the NLCG
        """
        super(NLCG, self).restart()
        self.NLCG.restart()



