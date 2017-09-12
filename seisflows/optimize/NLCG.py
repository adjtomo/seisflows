
import sys
import numpy as np

from seisflows.config import custom_import, ParameterError
from seisflows.plugins import optimize

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class NLCG(custom_import('optimize', 'base')):
    """ Nonlinear conjugate gradient method
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # line search algorithm
        if 'LINESEARCH' not in PAR:
            setattr(PAR, 'LINESEARCH', 'Bracket')

        # NLCG periodic restart interval
        if 'NLCGMAX' not in PAR:
            setattr(PAR, 'NLCGMAX', np.inf)

        # NLCG conjugacy restart threshold
        if 'NLCGTHRESH' not in PAR:
            setattr(PAR, 'NLCGTHRESH', np.inf)

        super(NLCG, self).check()


    def setup(self):
        super(NLCG, self).setup()

        self.NLCG = getattr(optimize, 'NLCG')(
            path=PATH.OPTIMIZE,
            maxiter=PAR.NLCGMAX,
            thresh=PAR.NLCGTHRESH,
            precond=self.precond)


    def compute_direction(self):
        g_new = self.load('g_new')
        p_new, self.restarted = self.NLCG()
        self.save('p_new', p_new)


    def restart(self):
        super(NLCG, self).restart()
        self.NLCG.restart()



