
import os
import numpy as np

from seisflows.tools import unix

from seisflows.tools.math import dot
from seisflows.tools.tools import loadtxt, savetxt, loadnpy, savenpy



class NLCG:
    """ Nonlinear conjugate gradient method
    """

    def __init__(self, path='.', load=loadnpy, save=savenpy, thresh=1., maxiter=np.inf, precond=None):
        self.path = path
        self.load = load
        self.save = save
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
        g_new = self.load('g_new')

        if self.iter == 1:
            return -g_new, 0

        elif self.iter > self.maxiter:
            print 'restarting NLCG... [periodic restart]'
            self.restart()
            return -g_new, 1

        # compute search direction
        g_old = self.load('g_old')
        p_old = self.load('p_old')

        if self.precond:
            beta = pollak_ribere(g_new, g_old, self.precond)
            p_new = -self.precond(g_new) + beta*p_old
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

def fletcher_reeves(g_new, g_old, precond=lambda x : x):
    num = dot(precond(g_new), g_new)
    den = dot(g_old, g_old)
    beta = num/den
    return beta

def pollak_ribere(g_new, g_old, precond=lambda x : x):
    num = dot(precond(g_new), g_new-g_old)
    den = dot(g_old, g_old)
    beta = num/den
    return beta

def check_conjugacy(g_new, g_old):
    return abs(dot(g_new, g_old) / dot(g_new, g_new))

def check_descent(p_new, g_new):
    return dot(p_new, g_new) / dot(g_new, g_new)




