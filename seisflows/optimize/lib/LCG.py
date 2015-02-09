
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt

from seisflows.optimize.lib.LBFGS import LBFGS


class LCG:
    """ Linear conjugate gradient method
    """
    def __init__(self, path, thresh, lcgmax, precond_type=0):
        self.path = path
        unix.mkdir(self.path+'/'+'LCG')

        self.ilcg = 0
        self.iter = 0
        self.thresh = thresh
        self.lcgmax = lcgmax

        self.precond_type = precond_type
        if precond_type in [1, 2]:
            self.LBFGS = LBFGS(path, 3)


    def precond(self, r):
        if self.precond_type == 0:
            return r

        elif self.precond_type == 1:
            if self.iter == 1:
                #print ' made it here-1'
                return r

            elif self.ilcg == 0:
                #print ' made it here-2'
                self.LBFGS.update()

            y = self.LBFGS.solve(r)
            return y


    def initialize(self):
        unix.cd(self.path)
        self.iter += 1
        self.ilcg = 0

        r = loadnpy('g_new')
        x = np.zeros(r.size)
        y = self.precond(r)
        p = -y

        unix.cd(self.path)
        savenpy('LCG/x', x)
        savenpy('LCG/r', r)
        savenpy('LCG/y', y)
        savenpy('LCG/p', p)
        savetxt('LCG/ry', np.dot(r, y))

        return p


    def update(self, ap):
        unix.cd(self.path)
        self.ilcg += 1

        x = loadnpy('LCG/x')
        r = loadnpy('LCG/r')
        y = loadnpy('LCG/y')
        p = loadnpy('LCG/p')
        ry = loadtxt('LCG/ry')

        pap = np.dot(p, ap)
        if pap < 0:
            print ' Newton failed [negative curvature]'
            isdone = True
            return p, isdone
                       

        alpha = ry/pap
        x = x + alpha*p
        r = r + alpha*ap

        # check status
        if self.ilcg == self.lcgmax:
            isdone = True
        elif np.linalg.norm(r) > self.thresh:
            isdone = True
        else:
            isdone = False
        if isdone:
            return x, isdone

        # apply preconditioner
        y = self.precond(r)

        ry_old = ry
        ry = np.dot(r, y)
        beta = ry/ry_old
        p = -y + beta*p

        savenpy('LCG/x', x)
        savenpy('LCG/r', r)
        savenpy('LCG/y', y)
        savenpy('LCG/p', p)
        savetxt('LCG/ry', np.dot(r, y))

        return p, isdone
