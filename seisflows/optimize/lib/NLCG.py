
import os
import numpy as np

from seisflows.tools import unix

from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt



class NLCG:
    """ Nonlinear conjugate gradient method
    """

    def __init__(self, path='.', thresh=1., maxiter=np.inf, precond=False):
        self.path = path
        self.maxiter = maxiter
        self.thresh = thresh
        self.precond = precond

        try:
            self.iter = loadtxt(self.path+'/'+'NLCG/iter')
        except IOError:
            unix.mkdir(self.path+'/'+'NLCG')
            self.iter = 0


    def __call__(self):
        """ Returns NLCG search direction
        """
        self.iter += 1
        savetxt(self.path+'/'+'NLCG/iter', self.iter)

        unix.cd(self.path)
        g_new = loadnpy('g_new')

        if self.iter == 1:
            return -g_new, 0

        elif self.iter > self.maxiter:
            print 'restarting NLCG... [periodic restart]'
            self.restart()
            return -g_new, 1

        # compute search direction
        g_old = loadnpy('g_old')
        p_old = loadnpy('p_old')

        if self.precond:
            d = loadnpy('precond')
            beta = pollak_ribere(d**0.5 * g_new, d**0.5 * g_old)
            p_new = -d*g_new + beta*p_old
        else:
            beta = pollak_ribere(g_new, g_old)
            p_new = -g_new + beta*p_old

        # check restart conditions
        if check_conjugacy(g_new, g_old) > self.thresh:
            print 'restarting NLCG... [loss of conjugacy]'
            self.restart()
            return -g_new, 1

        elif check_descent(p_new, g_new) > 0.:
            print 'restarting NLCG... [not a descent direction]'
            self.restart()
            return -g_new, 1

        else:
            return p_new, 0


    def restart(self):
        """ Restarts algorithm
        """
        self.iter = 1
        savetxt(self.path+'/'+'NLCG/iter', self.iter)



### utility functions

def fletcher_reeves(g_new, g_old):
    num = np.dot(g_new, g_new)
    den = np.dot(g_old, g_old)
    beta = num/den
    return beta

def pollak_ribere(g_new, g_old):
    num = np.dot(g_new, g_new-g_old)
    den = np.dot(g_old, g_old)
    beta = num/den
    return beta

def check_conjugacy(g_new, g_old):
    return abs(np.dot(g_new, g_old) / np.dot(g_new, g_new))

def check_descent(p_new, g_new):
    return np.dot(p_new, g_new) / np.dot(g_new, g_new)




