
import sys
import numpy as np

from seisflows.config import custom_import, ParameterError
from seisflows.plugins import optimize

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class LBFGS(custom_import('optimize', 'base')):
    """ Limited memory BFGS algorithm
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # line search algorithm
        if 'LINESEARCH' not in PAR:
            setattr(PAR, 'LINESEARCH', 'Backtrack')

        # LBFGS memory
        if 'LBFGSMEM' not in PAR:
            setattr(PAR, 'LBFGSMEM', 3)

        # LBFGS periodic restart interval
        if 'LBFGSMAX' not in PAR:
            setattr(PAR, 'LBFGSMAX', np.inf)

        # LBFGS angle restart threshold
        if 'LBFGSTHRESH' not in PAR:
            setattr(PAR, 'LBFGSTHRESH', 0.)

        super(LBFGS, self).check()


    def setup(self):
        super(LBFGS, self).setup()

        self.LBFGS = getattr(optimize, 'LBFGS')(
            path=PATH.OPTIMIZE,
            memory=PAR.LBFGSMEM,
            maxiter=PAR.LBFGSMAX,
            thresh=PAR.LBFGSTHRESH,
            precond=self.precond)


    def compute_direction(self):
        g_new = self.load('g_new')
        p_new, self.restarted = self.LBFGS()
        self.save('p_new', p_new)


    def restart(self):
        super(LBFGS, self).restart()
        self.LBFGS.restart()

