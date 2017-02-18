
import sys
import numpy as np

from seisflows.config import custom_import, ParameterError
from seisflows.optimize.lib.NLCG import NLCG as lib

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class NLCG(custom_import('optimize', 'steepest_descent')):
    """ Nonlinear conjugate gradient method
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # NLCG memory
        if 'NLCGMEM' not in PAR:
            setattr(PAR, 'NLCGMEM', 3)

        # NLCG periodic restart interval
        if 'NLCGMAX' not in PAR:
            setattr(PAR, 'NLCGMAX', np.inf)

        # NLCG angle restart threshold
        if 'NLCGTHRESH' not in PAR:
            setattr(PAR, 'NLCGTHRESH', np.inf)

        super(NLCG, self).check()


    def setup(self):
        super(NLCG, self).setup()

        self.NLCG = lib(
            path=PATH.OPTIMIZE,
            maxiter=PAR.NLCGMAX,
            thresh=PAR.NLCGTHRESH,
            precond=self.precond())


    def compute_direction(self):
        g_new = self.load('g_new')
        p_new, self.restarted = self.NLCG()
        self.save('p_new', p_new)
        self.savetxt('s_new', self.dot(g_new, p_new))
        return p_new


    def restart(self):
        super(NLCG, self).restart()
        self.NLCG.restart()



