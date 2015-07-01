
import os
import pickle
import time

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import savetxt, exists


class LBFGS(object):
    """ Limited-memory BFGS algorithm

        Includes optional safeguards: periodic restarting and descent
        conditions.

        To conserve memory, most vectors are read from disk rather than 
        passed from a calling routine.
    """

    def __init__(self, path='.', mem=5, thresh=0., maxiter=np.inf):
        assert exists(path)
        unix.cd(path)
        unix.mkdir('LBFGS')

        self.path = path
        self.thresh = thresh
        self.iter = 0
        self.maxiter = maxiter

        self.kk = 0
        self.nn = mem


    def __call__(self):
        """ Returns L-BFGS search direction
        """
        self.iter += 1
        savetxt('LBFGS/iter', self.iter)

        g = loadnpy('g_new')
        if self.iter == 1:
            return -g, 1

        elif self.iter > self.maxiter:
            print 'restarting LBFGS... [periodic restart]'
            self.restart()
            return -g, 0

        S, Y = self.update()
        q = self.apply(g, S, Y)

        if self.check_status(g, q) != 0:
            print 'restarting LBFGS... [not a descent direction]'
            self.restart()
            return -g, 1

        else:
            return -q, 0


    def update(self):
        """ Updates L-BFGS algorithm history
        """
        unix.cd(self.path)

        self.kk += 1
        self.kk = min(self.kk, self.nn)

        s = loadnpy('m_new') - loadnpy('m_old')
        y = loadnpy('g_new') - loadnpy('g_old')
        d = len(s)

        if self.kk == 1:
            S = np.memmap('LBFGS/S', mode='w+', dtype='float32', shape=(d, self.nn))
            Y = np.memmap('LBFGS/Y', mode='w+', dtype='float32', shape=(d, self.nn))
            S[:, 0] = s
            Y[:, 0] = y
        else:
            S = np.memmap('LBFGS/S', mode='r+', dtype='float32', shape=(d, self.nn))
            Y = np.memmap('LBFGS/Y', mode='r+', dtype='float32', shape=(d, self.nn))
            S[:, 1:] = S[:, :-1]
            Y[:, 1:] = Y[:, :-1]
            S[:, 0] = s
            Y[:, 0] = y

        return S, Y


    def apply(self, q, S, Y):
        """ Applies L-BFGS inverse Hessian to given vector
        """
        unix.cd(self.path)

        # first matrix product
        kk = self.kk
        rh = np.zeros(kk)
        al = np.zeros(kk)
        for ii in range(0, kk):
            rh[ii] = 1/np.dot(Y[:,ii], S[:,ii])
            al[ii] = rh[ii]*np.dot(S[:,ii], q)
            q = q - al[ii]*Y[:,ii]

        # use scaled identity matrix to initialize inverse Hessian, as proposed
        # by Liu and Nocedal 1989
        sty = np.dot(Y[:,0], S[:,0])
        yty = np.dot(Y[:,0], Y[:,0])
        r = sty/yty*q

        # second matrix product
        for ii in range(kk-1, -1, -1):
            be = rh[ii]*np.dot(Y[:,ii], r)
            r = r + S[:,ii]*(al[ii] - be)

        return r


    def restart(self):
        """ Discards history and resets counters
        """
        self.iter = 0
        savetxt('LBFGS/iter', 0)

        unix.cd(self.path)
        S = np.memmap('LBFGS/S', mode='r+')
        Y = np.memmap('LBFGS/Y', mode='r+')
        S[:] = 0.
        Y[:] = 0.

        time.sleep(2)


    def check_status(self, g, r):
        theta = np.dot(g,r)/(np.dot(g,g)*np.dot(r,r))**0.5
        if theta < 0.:
            print 'restarting LBFGS... [not a descent direction]'
            return -1
        elif theta < self.thresh:
            print 'restarting LBFGS... [practical safeguard]'
            return -2
        else:
            return 0


