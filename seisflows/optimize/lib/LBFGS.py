
import os
import pickle
import time

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt, exists


class LBFGS(object):
    """ Limited-memory BFGS algorithm

        Includes optional safeguards: periodic restarting and descent
        conditions.

        To conserve memory, most vectors are read from disk rather than 
        passed from a calling routine.
    """

    def __init__(self, path='.', memory=5, itermax=np.inf, thresh=0.):
        self.path = path
        self.kmax = memory
        self.itermax = itermax
        self.thresh = thresh

        unix.cd(self.path)

        if exists('LBFGS/iter'):
            self.iter = int(loadtxt('LBFGS/iter'))

        else:
            self.iter = 0
            unix.mkdir('LBFGS')
            savetxt('LBFGS/iter', 0)


    def update(self):
        """ Updates L-BFGS inverse Hessian
        """
        unix.cd(self.path)

        k = min(self.iter-1, self.kmax)

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

        return S, Y


    def compute(self):
        """ Applies L-BFGS inverse Hessian to gradient
        """
        unix.cd(self.path)

        self.iter += 1
        savetxt('LBFGS/iter', self.iter)

        g = loadnpy('g_new')

        if self.iter == 1:
            return -g, 1

        elif self.iter > self.itermax:
            print 'restarting LBFGS... [periodic restart]'
            self.restart()
            return -g, 1

        # update algorithm history
        S, Y = self.update()

        q = loadnpy('g_new')
        k = min(self.iter-1, self.kmax)
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

        if self.check_status(g, r) != 0:
            print 'restarting LBFGS... [not a descent direction]'
            self.restart()
            return -g, 1

        else:
            self.restarted = False
            return -r, 0


    def restart(self):
        """ Discards history and resets counters
        """
        self.restarted = True
        self.iter = 0
        savetxt('LBFGS/iter', 0)

        unix.cd(self.path)
        S = np.memmap('LBFGS/S', mode='r+')
        Y = np.memmap('LBFGS/Y', mode='r+')
        S[:] = 0.
        Y[:] = 0.

        time.sleep(2)


    def check_status(self, g, r, require_descent=True):

        theta = np.dot(g,r)/(np.dot(g,g)*np.dot(r,r))**0.5
        #print ' theta: %e' % theta

        if not require_descent:
            return 0
        elif theta > self.thresh:
            return 0
        else:
            return -1


