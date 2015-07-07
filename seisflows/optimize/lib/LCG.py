
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt

from seisflows.optimize.lib.LBFGS import LBFGS


class LCG(object):
    """ CG solver
    """
    def __init__(self, path, thresh=np.inf, maxiter=np.inf, precond=None):
        self.path = path
        self.maxiter = maxiter
        self.precond = precond

        self.ilcg = 0
        self.iter = 0


    def initialize(self):
        unix.mkdir(self.path+'/'+'LCG')
        unix.cd(self.path)

        self.iter += 1
        self.ilcg = 0

        r = loadnpy('g_new')
        x = np.zeros(r.size)
        savenpy('LCG/x', x)
        savenpy('LCG/r', r)

        y = self.apply_precond(r)
        p = -y
        savenpy('LCG/y', y)
        savenpy('LCG/p', p)
        savetxt('LCG/ry', np.dot(r, y))


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
            print ' Stopping LCG [negative curvature]'
            isdone = True
            return isdone
                       
        alpha = ry/pap
        x += alpha*p
        r += alpha*ap
        savenpy('LCG/x', x)
        savenpy('LCG/r', r)

        # check status
        if self.check_status(ap) == 0:
            isdone = True
        elif self.ilcg >= self.maxiter:
            isdone = True
        else:
            isdone = False

        if not isdone:
            y = self.apply_precond(r)
            ry_old = ry
            ry = np.dot(r, y)
            beta = ry/ry_old
            p = -y + beta*p

            savenpy('LCG/y', y)
            savenpy('LCG/p', p)
            savetxt('LCG/ry', np.dot(r, y))

        return isdone


    ### dummy methods, can be overloaded

    def check_status(self, *args, **kwargs):
        return -1

    def apply_precond(self, r):
        return r


