
import os
import pickle
import time

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt, exists


class LBFGS:
    """ Limited memory BFGS
    """

    def __init__(self, path='.', kmax=5, itermax=np.inf, thresh=np.inf):
        self.path = path
        self.kmax = kmax
        self.itermax = itermax

        if itermax == 0:
            self.itermax = np.inf

        unix.cd(self.path)

        if exists('LBFGS/iter'):
            self.iter = int(loadtxt('LBFGS/iter'))
        else:
            self.iter = 0
            unix.mkdir('LBFGS')
            savetxt('LBFGS/iter', 0)


    def update(self):
        unix.cd(self.path)

        self.iter += 1
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

        savetxt('LBFGS/iter', self.iter)
        del S
        del Y


    def solve(self, path='g_new'):
        unix.cd(self.path)
        g = loadnpy('g_new')

        if self.iter > self.itermax:
            print 'restarting LBFGS... [perdioc restart]'
            self.restart()
            return -g

        # compute L-BFSGS direction
        q = loadnpy(path)
        n = len(q)
        S = np.memmap('LBFGS/S', mode='r', dtype='float32', shape=(n, self.kmax))
        Y = np.memmap('LBFGS/Y', mode='r', dtype='float32', shape=(n, self.kmax))
        k = min(self.iter, self.kmax)
        rh = np.zeros(k)
        al = np.zeros(k)

        for i in range(0, k):
            rh[i] = 1/np.dot(Y[:, i], S[:, i])
            al[i] = rh[i]*np.dot(S[:, i], q)
            q = q - al[i]*Y[:, i]

        # initialize inverse Hessian using scaling proposed by 
        # Liu and Nocedal 1989
        sty = np.dot(Y[:,0], S[:,0])
        yty = np.dot(Y[:,0], Y[:,0])
        r = sty/yty*q

        for i in range(k-1, -1, -1):
            be = rh[i]*np.dot(Y[:, i], r)
            r = r + S[:, i]*(al[i] - be)

        # check restart conditions
        if np.dot(g, -r)/np.dot(g, g) >= 0:
            print 'restarting LBFGS... [not a descent direction]'
            self.restart()
            return -g

        else:
            self.restarted = False
            return -r


    def restart(self):
        self.restarted = True
        self.iter = 0
        savetxt('LBFGS/iter', 0)

        unix.cd(self.path)
        S = np.memmap('LBFGS/S', mode='r+')
        Y = np.memmap('LBFGS/Y', mode='r+')
        S[:] = 0.
        Y[:] = 0.

        time.sleep(2)



