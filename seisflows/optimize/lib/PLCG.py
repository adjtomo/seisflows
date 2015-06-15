
from os.path import join

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.io import OutputWriter

from seisflows.optimize.lib.LBFGS import LBFGS as LBFGS_base
from seisflows.optimize.lib.LCG import LCG


class LBFGS(LBFGS_base):
    """ Adapts L-BFGS functionality from nonlinear optimization to
     preconditioning
    """

    def update(self):
        """ Updates L-BFGS inverse Hessian
        """
        unix.cd(self.path)

        self.iter += 1
        savetxt('LBFGS/iter', self.iter)

        k = min(self.iter, self.kmax)

        s = loadnpy('m_new') - loadnpy('m_old')
        y = loadnpy('g_new') - loadnpy('g_old')
        n = len(s)

        if k == 1:
            S = np.memmap('LBFGS/S', mode='w+', dtype='float32', shape=(n, self.kmax))
            Y = np.memmap('LBFGS/Y', mode='w+', dtype='float32', shape=(n, self.kmax))
            S[:, 0] = s
            Y[:, 0] = y

        else:
            S = np.memmap('LBFGS/S', mode='r+', dtype='float32', shape=(n, self.kmax))
            Y = np.memmap('LBFGS/Y', mode='r+', dtype='float32', shape=(n, self.kmax))
            S[:, 1:] = S[:, :-1]
            Y[:, 1:] = Y[:, :-1]
            S[:, 0] = s
            Y[:, 0] = y

        del S
        del Y


    def solve(self, path, require_descent=True):
        """ Applies L-BFGS inverse Hessian to vector read from given path
        """
        unix.cd(self.path)

        g = loadnpy('g_new')

        if self.iter > self.itermax:
            print 'restarting LBFGS... [periodic restart]'
            self.restart()
            return -g

        q = loadnpy(path)
        n = len(q)

        # load algorithm history
        S = np.memmap('LBFGS/S', mode='r', dtype='float32', shape=(n, self.kmax))
        Y = np.memmap('LBFGS/Y', mode='r', dtype='float32', shape=(n, self.kmax))

        k = min(self.iter, self.kmax)
        rh = np.zeros(k)
        al = np.zeros(k)

        # first matrix product
        for i in range(0, k):
            rh[i] = 1/np.dot(Y[:, i], S[:, i])
            al[i] = rh[i]*np.dot(S[:, i], q)
            q = q - al[i]*Y[:, i]

        # use scaled identity matrix to initialize inverse Hessian, as proposed
        # by Liu and Nocedal 1989
        sty = np.dot(Y[:,0], S[:,0])
        yty = np.dot(Y[:,0], Y[:,0])
        r = sty/yty*q

        # second matrix product
        for i in range(k-1, -1, -1):
            be = rh[i]*np.dot(Y[:, i], r)
            r = r + S[:, i]*(al[i] - be)
        del S
        del Y

        if self.check_status(g, r, require_descent) != 0:
            print 'restarting LBFGS... [not a descent direction]'
            self.restart()
            return -g

        else:
            self.restarted = False
            return -r


class PLCG(LCG):
    """ Preconditioned truncated-Newton CG solver

      Adds preconditioning and adaptive stopping to LCG base class
    """
    def __init__(self, path, eta=1., **kwargs):
        self.eta = eta

        super(PLCG, self).__init__(path, **kwargs)

        # prepare output writer
        self.logpath = join(path, 'output.plcg')
        self.writer = OutputWriter(self.logpath, width=14)


    def precond(self, r):
        if self.precond_type in ['LBFGS_3']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, memory=3)
                y = r
            elif self.ilcg == 0:
                self.LBFGS.update()
                y = -self.LBFGS.solve('LCG/r', require_descent=True)
            else:
                y = -self.LBFGS.solve('LCG/r', require_descent=False)
            return y

        elif self.precond_type in ['LBFGS_6']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, memory=6)
                y = r
            elif self.ilcg == 0:
                self.LBFGS.update()
                y = -self.LBFGS.solve('LCG/r', require_descent=True)
            else:
                y = -self.LBFGS.solve('LCG/r', require_descent=False)
            return y

        elif self.precond_type in ['LBFGS_9']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, memory=9)
                y = r
            elif self.ilcg == 0:
                self.LBFGS.update()
                y = -self.LBFGS.solve('LCG/r', require_descent=True)
            else:
                y = -self.LBFGS.solve('LCG/r', require_descent=False)
            return y

        elif self.precond_type in ['LBFGS_3_norestart']:
            if self.iter == 1:
                self.LBFGS = LBFGS(self.path, memory=3)
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
        r = loadnpy('LCG/r')

        LHS = _norm(g)
        RHS = _norm(r)

        # for comparison, calculates forcing term proposed by 
        # Eisenstat & Walker 1996
        try:
            g_new = _norm(g)
            g_old = _norm(loadnpy('g_old'))
            eta1996 = g_new/g_old
        except:
            eta1996 = 1.

        if verbose:
            print ' RATIO:', RHS/LHS
            print ''

        self.writer(
            self.iter,
            self.ilcg,
            LHS,
            RHS,
            LHS/RHS,
            eta1996)

        # check termination condition
        if LHS > eta * RHS:
            return 0
        else:
            return -1



### utility functions

def _norm(v):
    return float(np.linalg.norm(v))
