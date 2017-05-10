
import numpy as np

from seisflows.tools import unix
from seisflows.tools.tools import savetxt, exists, loadnpy, savenpy
from seisflows.tools.math import angle


class LBFGS(object):
    """ Limited-memory BFGS algorithm

        Includes optional safeguards: periodic restarting and descent
        conditions.

        To conserve memory, most vectors are read from disk rather than 
        passed from a calling routine.
    """

    def __init__(self, path='.', load=loadnpy, save=savenpy, memory=5, thresh=0., maxiter=np.inf, precond=None):
        assert exists(path)
        unix.cd(path)
        unix.mkdir('LBFGS')

        self.path = path
        self.load = load
        self.save = save
        self.thresh = thresh
        self.maxiter = maxiter
        self.precond = precond
        self.memory = memory

        self.iter = 0
        self.memory_used = 0


    def __call__(self):
        """ Returns L-BFGS search direction
        """
        self.iter += 1

        unix.cd(self.path)
        g = self.load('g_new')
        if self.iter == 1:
            return -g, 0

        elif self.iter > self.maxiter:
            print 'restarting LBFGS... [periodic restart]'
            self.restart()
            return -g, 1

        S, Y = self.update()
        q = self.apply(g, S, Y)

        status = self.check_status(g,q)
        if status != 0:
            self.restart()
            return -g, status
        else:
            return -q, status


    def update(self):
        """ Updates L-BFGS algorithm history
        """
        unix.cd(self.path)

        s = self.load('m_new') - self.load('m_old')
        y = self.load('g_new') - self.load('g_old')

        m = len(s)
        n = self.memory

        if self.memory_used == 0:
            S = np.memmap('LBFGS/S', mode='w+', dtype='float32', shape=(m, n))
            Y = np.memmap('LBFGS/Y', mode='w+', dtype='float32', shape=(m, n))
            S[:, 0] = s
            Y[:, 0] = y
            self.memory_used = 1

        else:
            S = np.memmap('LBFGS/S', mode='r+', dtype='float32', shape=(m, n))
            Y = np.memmap('LBFGS/Y', mode='r+', dtype='float32', shape=(m, n))
            S[:, 1:] = S[:, :-1]
            Y[:, 1:] = Y[:, :-1]
            S[:, 0] = s
            Y[:, 0] = y

            if self.memory_used < self.memory:
                self.memory_used += 1

        return S, Y


    def apply(self, q, S=[], Y=[]):
        """ Applies L-BFGS inverse Hessian to given vector
        """
        unix.cd(self.path)

        if S==[] or Y==[]:
            m = len(q)
            n = self.memory
            S = np.memmap('LBFGS/S', mode='w+', dtype='float32', shape=(m, n))
            Y = np.memmap('LBFGS/Y', mode='w+', dtype='float32', shape=(m, n))

        # first matrix product
        kk = self.memory_used
        rh = np.zeros(kk)
        al = np.zeros(kk)
        for ii in range(kk):
            rh[ii] = 1/np.dot(Y[:,ii], S[:,ii])
            al[ii] = rh[ii]*np.dot(S[:,ii], q)
            q = q - al[ii]*Y[:,ii]

        if self.precond:
            r = self.precond(q)
        else:
            r = q

        # use scaling M3 proposed by Liu and Nocedal 1989
        sty = np.dot(Y[:,0], S[:,0])
        yty = np.dot(Y[:,0], Y[:,0])
        r *= sty/yty

        # second matrix product
        for ii in range(kk-1, -1, -1):
            be = rh[ii]*np.dot(Y[:,ii], r)
            r = r + S[:,ii]*(al[ii] - be)

        return r


    def restart(self):
        """ Discards history and resets counters
        """
        self.iter = 1
        self.memory_used = 0

        unix.cd(self.path)
        S = np.memmap('LBFGS/S', mode='r+')
        Y = np.memmap('LBFGS/Y', mode='r+')
        S[:] = 0.
        Y[:] = 0.


    def check_status(self, g, r):
        theta = 180.*np.pi**-1*angle(g,r)
        if not 0. < theta < 90.:
            print 'restarting LBFGS... [not a descent direction]'
            return 1
        elif theta > 90. - self.thresh:
            print 'restarting LBFGS... [practical safeguard]'
            return 1
        else:
            return 0


