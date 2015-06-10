
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt

from seisflows.optimize.lib.LBFGS import LBFGS


class PLCG:
    """ Truncated-Newton CG solver

      Adds preconditioning and adaptive stopping conditions to LCG base class
    """
    def __init__(self, *args, **kwargs):
        super(self, PLCG).__init__.(*args, **kwargs)

        self.logpath = self.path +'/'+ 'output.plcg'

        self.writer = OutputWriter(self.logpath,
            ['nonlinear',
            'linear',
            'eta'])
            


    def precond(self, r):
        if self.precond_type in ['LBFGS_3']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, step_memory=3)
                y = r
            elif self.ilcg == 0:
                self.LBFGS.update()
                y = -self.LBFGS.solve('LCG/r', require_descent=True)
            else:
                y = -self.LBFGS.solve('LCG/r', require_descent=False)
            return y

        elif self.precond_type in ['LBFGS_6']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, step_memory=6)
                y = r
            elif self.ilcg == 0:
                self.LBFGS.update()
                y = -self.LBFGS.solve('LCG/r', require_descent=True)
            else:
                y = -self.LBFGS.solve('LCG/r', require_descent=False)
            return y

        elif self.precond_type in ['LBFGS_9']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, step_memory=9)
                y = r
            elif self.ilcg == 0:
                self.LBFGS.update()
                y = -self.LBFGS.solve('LCG/r', require_descent=True)
            else:
                y = -self.LBFGS.solve('LCG/r', require_descent=False)
            return y

        elif self.precond_type in ['LBFGS_3_norestart']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, step_memory=3)
                y = r
            elif self.ilcg == 0:
                self.LBFGS.update()
                y = -self.LBFGS.solve('LCG/r', require_descent=False)
            else:
                y = -self.LBFGS.solve('LCG/r', require_descent=False)
            return y

        else:
            return r


    def check_status(self, ap, eta=0.9999, verbose=True):
        """ Checks Eisenstat-Walker termination status
        """
        g = loadnpy('g_new')
        LHS = np.linalg.norm(g+ap)
        RHS = np.linalg.norm(g)

        # for comparison, calculate forcing term proposed by 
        # Eisenstat & Walker 1996
        if self.iter > 1:
            g_new = np.linalg.norm(g)
            g_old = np.linalg.norm(loadnpy('g_old'))
            g_ratio = g_new/g_old
        else:
            g_ratio = np.nan

        if verbose:
            # print numerical statistics
            print ' k+1/k:', g_ratio
            print ' LHS:  ', LHS
            print ' RHS:  ', RHS
            print ' RATIO:', LHS/RHS
            print ''

        # check termination condition
        if LHS < eta * RHS:
            return 0
        else:
            return 1

