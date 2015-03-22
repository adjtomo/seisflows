
import os
import pickle
import time

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt


class LBFGS:
    """ Limited memory BFGS
    """

    def __init__(self, path='.', kmax=5, iter=1):
        self.path = path
        self.kmax = kmax

        unix.cd(self.path)
        unix.mkdir('LBFGS')

        if iter == 1:
            savetxt('LBFGS/k', 0)


    def update(self):
        unix.cd(self.path)
        s = loadnpy('m_new') - loadnpy('m_old')
        y = loadnpy('g_new') - loadnpy('g_old')
        n = len(s)

        k = int(loadtxt('LBFGS/k'))
        k = min(k+1, self.kmax)

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

        savetxt('LBFGS/k', k)
        del S
        del Y


    def solve(self, path='g_new'):
        unix.cd(self.path)
        q = loadnpy(path)
        n = len(q)

        S = np.memmap('LBFGS/S', mode='r', dtype='float32', shape=(n, self.kmax))
        Y = np.memmap('LBFGS/Y', mode='r', dtype='float32', shape=(n, self.kmax))
        k = int(loadtxt('LBFGS/k'))

        rh = np.zeros(k)
        al = np.zeros(k)

        for i in range(0, k):
            rh[i] = 1/np.dot(Y[:, i], S[:, i])
            al[i] = rh[i]*np.dot(S[:, i], q)
            q = q - al[i]*Y[:, i]

        sty = np.dot(Y[:, 0], S[:, 0])
        yty = np.dot(Y[:, 0], Y[:, 0])
        r = sty/yty*q

        for i in range(k - 1, -1, -1):
            be = rh[i]*np.dot(Y[:, i], r)
            r = r + S[:, i]*(al[i] - be)

        # check for ill conditioning
        g = loadnpy('g_new')
        if np.dot(g, -r)/np.dot(g, g) >= 0:
            self.restart()
            r = g
            self.restarted = True
        else:
            self.restarted = False

        return -r


    def restart(self):
        print 'restarting LBFGS...'
        time.sleep(2)

        unix.cd(self.path)
        savetxt('LBFGS/k', 0)
        S = np.memmap('LBFGS/S', mode='r+')
        Y = np.memmap('LBFGS/Y', mode='r+')
        S[:] = 0.
        Y[:] = 0.


