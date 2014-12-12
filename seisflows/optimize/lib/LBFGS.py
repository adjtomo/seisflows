
import os
import pickle
import time

import numpy as np

from seisflows.tools import unix


class LBFGS:
    """ Limited memory BFGS
    """

    def __init__(self,path='.',kmax=5,iter=1):

        self.path = path
        self.load = load
        self.save = save
        self.kmax = kmax

        unix.mkdir(self.path+'/'+'LBFGS')
        unix.cd(self.path+'/'+'LBFGS')

        if iter == 1:
            savetxt('k',0)


    def update(self):

        unix.cd(self.path)
        s = self.load('m_new') - self.load('m_old')
        y = self.load('g_new') - self.load('g_old')
        n = len(s)

        unix.cd('LBFGS')
        k = loadtxt('k')
        k = min(k+1,self.kmax)

        if k == 1:
            S = np.memmap('S',mode='w+',dtype='float32',shape=(n,self.kmax))
            Y = np.memmap('Y',mode='w+',dtype='float32',shape=(n,self.kmax))
            S[:,0] = s
            Y[:,0] = y

        else:
            S = np.memmap('S',mode='r+',dtype='float32',shape=(n,self.kmax))
            Y = np.memmap('Y',mode='r+',dtype='float32',shape=(n,self.kmax))
            S[:,1:] = S[:,:-1]
            Y[:,1:] = Y[:,:-1]
            S[:,0] = s
            Y[:,0] = y

        savetxt('k',k)
        del S
        del Y


    def solve(self):

        unix.cd(self.path)
        g = self.load('g_new')
        q = np.copy(g)
        n = len(q)

        unix.cd('LBFGS')
        S = np.memmap('S',mode='r',dtype='float32',shape=(n,self.kmax))
        Y = np.memmap('Y',mode='r',dtype='float32',shape=(n,self.kmax))
        k = loadtxt('k')

        rh = np.zeros(k)
        al = np.zeros(k)

        for i in range(0,k):
            rh[i] = 1/np.dot(Y[:,i],S[:,i])
            al[i] = rh[i]*np.dot(S[:,i],q)
            q = q - al[i]*Y[:,i]

        sty = np.dot(Y[:,0],S[:,0])
        yty = np.dot(Y[:,0],Y[:,0])
        r = sty/yty * q

        for i in range(k-1,-1,-1):

            be = rh[i]*np.dot(Y[:,i],r)
            r = r + S[:,i]*(al[i]-be)

        # check for ill conditioning
        if np.dot(g,-r) >= 0:
            self.restart()
            return g

        return r


    def restart(self):

        print 'restarting LBFGS...'
        time.sleep(2)

        unix.cd(self.path+'/'+'LBFGS')
        savetxt('k',0)
        S = np.memmap('S',mode='r+')
        Y = np.memmap('Y',mode='r+')
        S[:] = 0.
        Y[:] = 0.



# utility functions

def loadtxt(filename):
    return int(np.loadtxt(filename))

def savetxt(filename,v):
    np.savetxt(filename,[v],'%d')

def load(filename):
    return np.load(filename)

def save(filename,v):
    np.save(filename,v)
    unix.mv(filename+'.npy',filename)
