
import numpy as np

from seisflows.tools import unix

import LBFGS


class LCG:
    """ Linear conjugate gradient method
    """

    def __init__(self, path, thresh, itermax, precond_type=0):
        self.path = path
        unix.mkdir(self.path+'/'+'LCG')

        self.load = load
        self.save = save

        self.iter = 0
        self.thresh = thresh
        self.itermax = itermax

        self.precond_type = precond_type
        if precond_type in [1, 2]:
            self.LBFGS = LBFGS(path+'/'+'LCG', load, save, itermax)

    def precond(self, r):
        # no preconditioner
        if self.precond_type == 0:
            y = r

        # LBFGS without restarts
        elif self.prcond_type == 1:
            if self.iter == 1:
                self.LBFGS = LBFGS()
                y = r
            else:
                self.LBFGS.update()
                y = self.LBFGS.solve(r)

        # LBFGS with restarts
        elif self.precond_type == 2:
            if self.iter == 1:
                y = r
            else:
                self.LBFGS.update()
                y = -self.LBFGS.solve()

        return y

    def initialize(self):

        unix.cd(self.path)
        self.iter = 0

        r = self.load('g_new')
        x = np.zeros(r.size)
        y = self.precond(r)
        p = -y

        self.save('LCG/x', x)
        self.save('LCG/r', r)
        self.save('LCG/y', y)
        self.save('LCG/p', p)
        savetxt('LCG/ry', np.dot(r, y))

        return p

    def update(self, ap):

        unix.cd(self.path)
        self.iter += 1

        x = self.load('LCG/x')
        r = self.load('LCG/r')
        y = self.load('LCG/y')
        p = self.load('LCG/p')
        ry = loadtxt('LCG/ry')

        pap = np.dot(p, ap)
        alpha = ry/pap
        x = x + alpha*p
        r = r + alpha*ap

        # check status
        if self.iter == self.itermax:
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

        self.save('LCG/x', x)
        self.save('LCG/r', r)
        self.save('LCG/y', y)
        self.save('LCG/p', p)
        savetxt('LCG/ry', np.dot(r, y))

        return p, isdone


# -- utility functions

def loadtxt(filename):
    return np.loadtxt(filename)


def savetxt(filename, v):
    np.savetxt(filename, [v], '%e')


def load(filename):
    return np.load(filename)


def save(filename, v):
    np.save(filename, v)
    unix.mv(filename+'.npy', filename)
