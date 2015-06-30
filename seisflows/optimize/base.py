
from os.path import join

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError
from seisflows.tools.io import OutputWriter

from seisflows.tools.math import polyfit2, backtrack2
from seisflows.optimize.lib.LBFGS import LBFGS
from seisflows.optimize.lib.NLCG import NLCG

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class base(object):
    """ Nonlinear optimization base class.

     Available nonlinear optimization algorithms include steepest descent (SD),
     nonlinear conjugate gradient (NLCG), and LBFGS. Associated step control 
     algorithms include a backtracking line search and a bracketing line search.

     While NLCG (a Krylov method) and LBFGS (a quasi-Newton metod) are both 
     widely used for geophysical inversion, LBFGS is more efficient and more
     robust. NLCG requires occasional restarts to avoid numerical stagnation. 
     LBFGS generally requires fewer restarts. Restarts are controlled by 
     numerical tuning parameters. Default values provided for these parameters 
     (see below) should work well for most geophyscial inversions.

     To reduce memory overhead, input vectors are read from the directory
     self.path rather than passed from a calling routine. At the start of each
     search direction computation, the current model and gradient are read from
     'm_new' and 'g_new'; the resulting search direction is written to 'p_new'.
     As the optimization procedure progresses, other information is stored in
     the self.path directory as well.
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        if 'BEGIN' not in PAR:
            raise ParameterError

        if 'END' not in PAR:
            raise ParameterError

        if 'SUBMIT' not in PATH:
            raise ParameterError

        if 'OPTIMIZE' not in PATH:
            setattr(PATH, 'OPTIMIZE', join(PATH.GLOBAL, 'optimize'))

        # search direction algorithm
        if 'SCHEME' not in PAR:
            setattr(PAR, 'SCHEME', 'LBFGS')

        # line search algorithm
        if 'LINESEARCH' not in PAR:
            if  PAR.SCHEME in ['LBFGS']:
                setattr(PAR, 'LINESEARCH', 'Backtrack')

            elif 'Newton' in PAR.SCHEME:
                setattr(PAR, 'LINESEARCH', 'Backtrack')

            else:
                setattr(PAR, 'LINESEARCH', 'Bracket')

        # search direction tuning parameters
        if 'NLCGMAX' not in PAR:
            setattr(PAR, 'NLCGMAX', np.inf)

        if 'NLCGTHRESH' not in PAR:
            setattr(PAR, 'NLCGTHRESH', np.inf)

        if 'LBFGSMEMORY' not in PAR:
            setattr(PAR, 'LBFGSMEMORY', 5)

        if 'LBFGSTHRESH' not in PAR:
            setattr(PAR, 'LBFGSTHRESH', 0.)

        # line search tuning paraemters
        if 'STEPMAX' not in PAR:
            setattr(PAR, 'STEPMAX', 10)

        if 'STEPTHRESH' not in PAR:
            setattr(PAR, 'STEPTHRESH', None)

        if 'STEPINIT' not in PAR:
            setattr(PAR, 'STEPINIT', 0.05)

        if 'STEPFACTOR' not in PAR:
            setattr(PAR, 'STEPFACTOR', 0.25)

        if 'STEPOVERSHOOT' not in PAR:
            setattr(PAR, 'STEPOVERSHOOT', 0.)


    def setup(self):
        """ Sets up nonlinear optimization machinery
        """
        # specify paths
        self.path = PATH.OPTIMIZE
        self.logpath = join(PATH.SUBMIT, 'output.optim')
        self.writer = OutputWriter(self.logpath)
        unix.mkdir(self.path)

        # prepare algorithm machinery
        if PAR.SCHEME in ['NLCG']:
            self.NLCG = NLCG(
                self.path,
                max=PAR.NLCGMAX,
                thresh=PAR.NLCGTHRESH)

        elif PAR.SCHEME in ['LBFGS']:
            self.LBFGS = LBFGS(
                self.path, 
                memory=PAR.LBFGSMEMORY, 
                thresh=PAR.LBFGSTHRESH)

        self.restart = 0
        self.restart_count = 0


    ### search direction methods

     # The following names are used in the search direction methods
     # and for writing data to self.path:
     #    m_new - current model
     #    m_old - previous model
     #    m_try - trial model
     #    f_new - current objective function value
     #    f_old - previous objective function value
     #    f_try - trial objective function value
     #    g_new - current gradient direction
     #    g_old - previous gradient direction
     #    p_new - current search direction
     #    p_old - previous search direction
     #    s_new - current slope along search direction
     #    s_old - previous slope along search direction
     #    alpha - trial step length


    def compute_direction(self):
        """ Computes model update direction from stored gradient
        """
        unix.cd(self.path)
        g_new = loadnpy('g_new')

        if PAR.SCHEME == 'SD':
            p_new = -g_new

        elif PAR.SCHEME == 'NLCG':
            p_new, self.restart = self.NLCG()

        elif PAR.SCHEME =='LBFGS':
            p_new, self.restart = self.LBFGS()

        # keep track of number of restarts
        if self.restart:
            self.restart_count += 1

        savenpy('p_new', p_new)
        savetxt('s_new', np.dot(g_new, p_new))


    ### line search methods

    # The following names are used in the line search methods:
    #     m - model vector
    #     p - search direction vector
    #     s - slope along search direction
    #     f - value of objective function, evaluated at m
    #     x - step length along search direction
    #     p_ratio - ratio of model norm to search direction norm
    #     s_ratio - ratio of current slope to previous slope

    def initialize_search(self):
        """ Determines initial step length for line search
        """
        unix.cd(self.path)

        m = loadnpy('m_new')
        p = loadnpy('p_new')
        f = loadtxt('f_new')
        norm_m = max(abs(m))
        norm_p = max(abs(p))
        p_ratio = float(norm_m/norm_p)

        # reset search history
        self.search_history = [[0., f]]
        self.step_count = 0
        self.isdone = 0
        self.isbest = 0
        self.isbrak = 0

        # determine initial step length
        if self.iter == 1:
            alpha = p_ratio * PAR.STEPINIT
        elif self.restart:
            alpha = self.step_init()
        elif PAR.LINESEARCH in ['Backtrack']:
            alpha = 1.
        else:
            alpha = self.step_init()

        # optional ad hoc scaling
        if PAR.STEPOVERSHOOT:
            alpha *= PAR.STEPOVERSHOOT

        # optional maximum step length safegaurd
        if PAR.STEPTHRESH:
            if alpha > p_ratio * PAR.STEPTHRESH:
                alpha = p_ratio * PAR.STEPTHRESH

        # prepare trial model
        savenpy('m_try', m + p*alpha)
        savetxt('alpha', alpha)

        # prepare output writer
        if self.iter == 1:
            self.writer.header('iter', 'steplen', 'misfit')
        self.writer(self.iter, 0., f)


    def search_status(self):
        """ Determines status of line search

            Maintains line search history by keeping track of step length and
            function value from each trial model evaluation. From line search
            history, determines whether stopping criteria have been satisfied.
        """
        unix.cd(self.path)

        # update search history
        x_ = loadtxt('alpha')
        f_ = loadtxt('f_try')

        if np.isnan(f_):
            raise ValueError

        self.search_history += [[x_, f_]]
        self.step_count += 1

        x = self.step_lens()
        f = self.func_vals()

        # is current step length the best so far?
        vals = self.func_vals(sort=False)
        if np.all(vals[-1] < vals[:-1]):
            self.isbest = 1

        # are stopping criteria satisfied?
        if PAR.LINESEARCH == 'Bracket' or\
           self.restart:
            if self.isbrak:
                self.isbest = 1
                self.isdone = 1
            elif any(f[1:] < f[0]) and (f[-2] < f[-1]):
                self.isbrak = 1

        elif PAR.LINESEARCH == 'Backtrack':
            if any(f[1:] < f[0]):
                self.isdone = 1

        elif PAR.LINESEARCH == 'Fixed':
            if any(f[1:] < f[0]) and (f[-2] < f[-1]):
                self.isdone = 1

        # append latest line search values
        self.writer(None, x_, f_)
        if self.isdone:
            self.writer.newline()

        if self.step_count >= PAR.STEPMAX:
            self.isdone = -1
            print ' line search failed [max iter]\n'

        return self.isdone, self.isbest


    def compute_step(self):
        """ Computes next trial step length
        """
        unix.cd(self.path)
        m = loadnpy('m_new')
        p = loadnpy('p_new')
        s = loadtxt('s_new')

        norm_m = max(abs(m))
        norm_p = max(abs(p))
        p_ratio = float(norm_m/norm_p)

        x = self.step_lens()
        f = self.func_vals()

        # compute trial step length
        if PAR.LINESEARCH == 'Bracket' or\
           self.restart:
            if any(f[1:] < f[0]) and (f[-2] < f[-1]):
                alpha = polyfit2(x, f)
            elif any(f[1:] < f[0]):
                alpha = loadtxt('alpha')*PAR.STEPFACTOR**-1
            else:
                alpha = loadtxt('alpha')*PAR.STEPFACTOR

        elif PAR.LINESEARCH == 'Backtrack':
            alpha = backtrack2(f[0], s, x[1], f[1], b1=0.1, b2=0.5)

        elif PAR.LINESEARCH == 'Fixed':
            alpha = p_ratio*(step + 1)*PAR.STEPINIT

        # write trial model
        savetxt('alpha', alpha)
        savenpy('m_try', m + p*alpha)


    def finalize_search(self):
        """ Cleans working directory and writes updated model
        """
        unix.cd(self.path)
        m = loadnpy('m_new')
        p = loadnpy('p_new')

        x = self.step_lens()
        f = self.func_vals()

        # clean working directory
        unix.rm('alpha')
        unix.rm('m_try')
        unix.rm('f_try')

        if self.iter > 1:
            unix.rm('m_old')
            unix.rm('f_old')
            unix.rm('g_old')
            unix.rm('p_old')
            unix.rm('s_old')

        unix.mv('m_new', 'm_old')
        unix.mv('f_new', 'f_old')
        unix.mv('g_new', 'g_old')
        unix.mv('p_new', 'p_old')
        unix.mv('s_new', 's_old')

        # write updated model
        alpha = x[f.argmin()]
        savetxt('alpha', alpha)
        savenpy('m_new', m + p*alpha)
        savetxt('f_new', f.min())


    ### line search utilities

    def step_init(self):
        alpha = loadtxt('alpha')
        s_new = loadtxt('s_new')
        s_old = loadtxt('s_old')
        s_ratio = s_new/s_old
        return 2.*s_ratio*alpha


    def step_lens(self, sort=True):
        x, f = zip(*self.search_history)
        x = np.array(x)
        f = np.array(f)
        f_sorted = f[abs(x).argsort()]
        x_sorted = x[abs(x).argsort()]
        if sort:
            return x_sorted
        else:
            return x


    def func_vals(self, sort=True):
        x, f = zip(*self.search_history)
        x = np.array(x)
        f = np.array(f)
        f_sorted = f[abs(x).argsort()]
        x_sorted = x[abs(x).argsort()]
        if sort:
            return f_sorted
        else:
            return f
