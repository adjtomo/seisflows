#!/usr/bin/env python
"""
This is the custom class for an LBFGS optimization schema.
It supercedes the `seisflows.optimize.base` class
"""
import sys
import numpy as np

from seisflows.optimize.base import Base
from seisflows.plugins import optimize

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class LBFGS(Base):
    """
    The Limited memory BFGS algorithm
    Calls upon seisflows.plugin.optimize.LBFGS to accomplish LBFGS algorithm
    """
    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings

        :type LBFGS: Class
        :param LBFGS: plugin LBFGS class that controls the machinery of the
            L-BFGS optimization schema
        """
        self.LBFGS = None

    def check(self):
        """
        Checks parameters, paths, and dependencies
        """
        # Line search algorithm
        if "LINESEARCH" not in PAR:
            setattr(PAR, "LINESEARCH", "Backtrack")

        # LBFGS memory
        if "LBFGSMEM" not in PAR:
            setattr(PAR, "LBFGSMEM", 3)

        # LBFGS periodic restart interval
        if "LBFGSMAX" not in PAR:
            setattr(PAR, "LBFGSMAX", np.inf)

        # LBFGS angle restart threshold
        if "LBFGSTHRESH" not in PAR:
            setattr(PAR, "LBFGSTHRESH", 0.)

        # Include all checks from Base class
        super(LBFGS, self).check()

    def setup(self):
        """
        Set up the LBFGS optimization schema
        """
        super(LBFGS, self).setup()
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
        super(LBFGS, self).restart()
        self.LBFGS.restart()

