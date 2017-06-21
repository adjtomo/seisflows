
from os.path import join

import numpy as np

from seisflows.tools import unix
from seisflows.tools.tools import loadtxt, savetxt
#from seisflows.tools.io import OutputWriter

from seisflows.plugins.optimize.LBFGS import LBFGS
from seisflows.plugins.optimize.LCG import LCG
 
 
class LBFGS_(LBFGS):
    """ Adapts L-BFGS from nonlinear optimization to preconditioning
    """
    pass



class PLCG(LCG):
    """ Preconditioned truncated-Newton CG solver

      Adds preconditioning and adaptive stopping to LCG base class
    """
    def __init__(self, path, eta=1., **kwargs):
        self.eta = eta

        super(PLCG, self).__init__(path, **kwargs)

        # prepare output writer
        self.logpath = join(path, 'output.plcg')
        #self.writer = OutputWriter(self.logpath, width=14)


    def apply_precond(self, r):
        if not self.precond:
            return r

        elif self.precond in ['LBFGS_3']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, memory=3)
                y = r
            elif self.ilcg == 0:
                S, Y = self.LBFGS.update()
                y = -self.LBFGS.apply(self.load('LCG/r'), S, Y)
            else:
                y = -self.LBFGS.apply(self.load('LCG/r'))
            return y

        elif self.precond in ['LBFGS_6']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, memory=6)
                y = r
            elif self.ilcg == 0:
                S, Y = self.LBFGS.update()
                y = -self.LBFGS.apply(self.load('LCG/r'), S, Y)
            else:
                y = -self.LBFGS.apply(self.load('LCG/r'))
            return y

        elif self.precond in ['LBFGS_9']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, memory=9)
                y = r
            elif self.ilcg == 0:
                S, Y = self.LBFGS.update()
                y = -self.LBFGS.apply(self.load('LCG/r'), S, Y)
            else:
                y = -self.LBFGS.apply(self.load('LCG/r'))
            return y

        else:
            raise ValueError



    def check_status(self, ap, verbose=True):
        """ Checks Eisenstat-Walker termination status
        """
        g0 = self.load('g_new')
        g1 = self.load('LCG/r')

        LHS = _norm(g1)
        RHS = _norm(g0)

        # for comparison, calculates forcing term proposed by 
        # Eisenstat & Walker 1996
        try:
            g_new = _norm(g)
            g_old = _norm(self.load('g_old'))
            eta1996 = g_new/g_old
        except:
            eta1996 = 1.

        if verbose:
            print ' RATIO:', LHS/RHS
            print ''

        #self.writer(
        #    self.iter,
        #    self.ilcg,
        #    LHS,
        #    RHS,
        #    LHS/RHS,
        #    eta1996)

        # check termination condition
        if LHS < self.eta * RHS:
            return _done
        else:
            return not _done



### utility functions

_done = 0

def _norm(v):
    return float(np.linalg.norm(v))
