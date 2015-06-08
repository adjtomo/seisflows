
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt

from seisflows.optimize.lib.LBFGS import LBFGS


class LCG:
    """ Truncated-Newton CG solver
    """
    def __init__(self, path, eta0, lcgmax, precond_type=None):
        self.path = path
        unix.mkdir(self.path+'/'+'LCG')

        self.ilcg = 0
        self.iter = 0
        self.eta0 = eta0
        self.lcgmax = lcgmax
        self.precond_type = precond_type


    def precond(self, r):
        if self.precond_type in ['NonlinearLBFGS']:
            if self.iter == 1:
                self.LBFGS = LBFGS(path, 3)
                y = r
            elif self.ilcg == 0:
                self.LBFGS.update()
                y = -self.LBFGS.solve('LCG/r')
            else:
                y = -self.LBFGS.solve('LCG/r')
            return y

        else:
            return r


    def initialize(self):
        unix.cd(self.path)
        self.iter += 1
        self.ilcg = 0

        r = loadnpy('g_new')
        x = np.zeros(r.size)
        savenpy('LCG/x', x)
        savenpy('LCG/r', r)

        y = self.precond(r)
        p = -y
        savenpy('LCG/y', y)
        savenpy('LCG/p', p)
        savetxt('LCG/ry', np.dot(r, y))


    def update(self, ap):
        """ performs CG update
        """

        # utility function
        def EisenstatWalker(eta=.999, verbose=True):

            g = loadnpy('g_new')
            if self.iter > 1:
                g_new = np.linalg.norm(g)
                g_old = np.linalg.norm(loadnpy('g_old'))
                g_ratio = g_new/g_old

            LHS = np.linalg.norm(g+ap)
            RHS = np.linalg.norm(g)

            if verbose:
                print ' LHS:  ', LHS
                print ' RHS:  ', RHS
                print ' RATIO:', LHS/RHS
                print ''

            if LHS < eta * RHS:
                isdone = True
            else:
                isdone = False
            return isdone


        unix.cd(self.path)
        self.ilcg += 1

        x = loadnpy('LCG/x')
        r = loadnpy('LCG/r')
        y = loadnpy('LCG/y')
        p = loadnpy('LCG/p')
        ry = loadtxt('LCG/ry')

        pap = np.dot(p, ap)
        if pap < 0:
            print ' Stopping LCG [negative curvature]'
            isdone = True
            return isdone
                       
        alpha = ry/pap
        x = x + alpha*p
        r = r + alpha*ap
        savenpy('LCG/x', x)
        savenpy('LCG/r', r)

        # check status
        if self.ilcg >= self.lcgmax:
            isdone = True
        elif EisenstatWalker():
            isdone = True
        else:
            isdone = False

        if not isdone:
            y = self.precond(r)
            ry_old = ry
            ry = np.dot(r, y)
            beta = ry/ry_old
            p = -y + beta*p

            savenpy('LCG/y', y)
            savenpy('LCG/p', p)
            savetxt('LCG/ry', np.dot(r, y))

        return isdone
