
import os
import numpy as np

from seisflows.tools import unix


class NLCG:
    """ Nonlinear conjugate gradient method
    """

    def __init__(self, path, thresh, itercgmax):

        self.path = path
        unix.mkdir(self.path+'/'+'NLCG')

        self.load = load
        self.save = save

        self.itercgmax = itercgmax
        self.thresh = thresh

        try:
            self.itercg = loadtxt(self.path+'/'+'NLCG/itercg')
        except IOError:
            self.itercg = 0

    def compute(self):

        # utility functions

        def fletcher_reeves():
            # fletcher-reeves update
            top = np.dot(g_new, g_new)
            bot = np.dot(g_old, g_old)
            beta = top/bot
            return -g_new + beta*p_old

        def pollak_ribere():
            # pollack-riebere update
            top = np.dot(g_new, g_new-g_old)
            bot = np.dot(g_old, g_old)
            beta = top/bot
            return -g_new + beta*p_old

        def rho():
            return abs(np.dot(g_new, g_old) / np.dot(g_new, g_new))


        # algorithm starts here
        self.itercg += 1
        savetxt(self.path+'/'+'NLCG/itercg', self.itercg)

        unix.cd(self.path)
        g_new = self.load('g_new')

        if self.itercg == 1:
            return -g_new

        if self.itercg > self.itercgmax:
            print 'restarting NLCG... [periodic restart]'
            self.restart()
            return -g_new

        # compute NLCG direction
        g_old = self.load('g_old')
        p_old = self.load('p_old')
        p_new = pollak_ribere()

        # check restart conditions
        if rho() > self.thresh:
            print 'restarting NLCG... [loss of conjugacy]'
            self.restart()
            return -g_new

        elif np.dot(p_new, g_new) > 0:
            print 'restarting NLCG... [not a descent direction]'
            self.restart()
            return -g_new

        else:
            return p_new


    def restart(self):
        self.itercg = 1
        savetxt(self.path+'/'+'NLCG/itercg', self.itercg)



def loadtxt(filename):
    return int(np.loadtxt(filename))


def savetxt(filename, v):
    np.savetxt(filename, [v], '%d')


def load(filename):
    return np.load(filename)


def save(filename, v):
    np.save(filename, v)
    unix.mv(filename+'.npy', filename)
