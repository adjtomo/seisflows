
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

        def fletcher_reeves():
            # fletcher-reeves update
            top = np.dot(g_new, g_new)
            bot = np.dot(g_old, g_old)
            beta = top/bot
            return -g_new + beta*p_old

        def pollak_ribere():
            # pollack-riebere update
            top = np.dot(g_new, g_new - g_old)
            bot = np.dot(g_old, g_old)
            beta = top/bot
            return -g_new + beta*p_old

        def rho():
            return abs(np.dot(g_new, g_old) / np.dot(g_new, g_new))

        self.itercg += 1
        unix.cd(self.path)

        if self.itercg == 1:
            g_new = self.load('g_new')

            unix.cd('NLCG')
            p_new = -g_new

        elif self.itercg > 1:
            g_new = self.load('g_new')
            g_old = self.load('g_old')
            p_old = self.load('p_old')

            unix.cd('NLCG')

            # FIXME: do the following two 'if' clauses need to be separated ?
            #        is it important to know that in some case we have periodic
            #        restart AND loss of conjugacy?
            if self.itercg > self.itercgmax:
                # require periodic restarts
                print 'restarting NLCG... [periodic restart]'
                self.itercg = 1
                p_new = -g_new

            else:
                p_new = pollak_ribere()

            if rho() > self.thresh:
                # require orthogonality
                print 'restarting NLCG... [loss of conjugacy]'
                self.itercg = 1
                p_new = -g_new

            elif np.dot(p_new, g_new) > 0:
                # require descent direction
                print 'restarting NLCG... [not a descent direction]'
                self.itercg = 1
                p_new = -g_new

        savetxt(self.path+'/'+'NLCG/itercg', self.itercg)
        return p_new


def loadtxt(filename):
    return int(np.loadtxt(filename))


def savetxt(filename, v):
    np.savetxt(filename, [v], '%d')


def load(filename):
    return np.load(filename)


def save(filename, v):
    np.save(filename, v)
    unix.mv(filename+'.npy', filename)
