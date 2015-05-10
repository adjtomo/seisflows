
from os.path import join

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError
from seisflows.tools.io import OutputWriter
from seisflows.optimize import lib

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class base(object):
    """ Nonlinear optimization base class.

     Available nonlinear optimization algorithms include gradient descent,
     nonlinear conjugate gradient method, and limited-memory BFGS 
     (a quasi-Newton method).

     Available line search algorithms include a backtracking line search based
     on quadratic interpolation and a bracketing and interpolation procedure
     (abbreviated as 'Backtrack' and 'Bracket' respectively.)

     To reduce memory overhead, input vectors are read from the directory
     cls.path rather than passed from a calling routine. At the start of each
     search direction computation, the current model and gradient are read from
     'm_new' and 'g_new'; the resulting search direction is written to 'p_new'.
     As the optimization procedure progresses, other information is stored in
     the cls.path directory as well.
    """

    def check(cls):
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

        # search direction parameters
        if 'SCHEME' not in PAR:
            setattr(PAR, 'SCHEME', 'QuasiNewton')

        if 'NLCGMAX' not in PAR:
            setattr(PAR, 'NLCGMAX', 10)

        if 'NLCGTHRESH' not in PAR:
            setattr(PAR, 'NLCGTHRESH', 1.0)

        if 'LBFGSMAX' not in PAR:
            setattr(PAR, 'LBFGSMAX', 6)

        # line search parameters
        if 'LINESEARCH' not in PAR:
            if 'Newton' in PAR.SCHEME:
                setattr(PAR, 'LINESEARCH', 'Backtrack')
            else:
                setattr(PAR, 'LINESEARCH', 'Bracket')

        if 'SRCHMAX' not in PAR:
            setattr(PAR, 'SRCHMAX', 10)

        if 'STEPLEN' not in PAR:
            setattr(PAR, 'STEPLEN', 0.05)

        if 'STEPLENMAX' not in PAR:
            setattr(PAR, 'STEPLENMAX', None)

        if 'ADHOCSCALING' not in PAR:
            setattr(PAR, 'ADHOCSCALING', 0.)



    def setup(cls):
        """ Sets up directory in which to store optimization vectors
        """
        cls.path = PATH.OPTIMIZE
        unix.mkdir(cls.path)

        # prepare algorithm machinery
        if PAR.SCHEME in ['ConjugateGradient']:
            cls.NLCG = lib.NLCG(cls.path, PAR.NLCGTHRESH, PAR.NLCGMAX)

        elif PAR.SCHEME in ['QuasiNewton']:
            cls.LBFGS = lib.LBFGS(cls.path, PAR.LBFGSMAX, PAR.BEGIN)

        # prepare output writer
        cls.writer = OutputWriter(PATH.SUBMIT + '/' + 'output.optim',
            ['iter', 'step', 'misfit'])


    ### search direction methods

     # The following names are used in the search direction methods
     # and for writing data to cls.path:
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


    def compute_direction(cls):
        """ Computes model update direction from stored gradient
        """
        unix.cd(cls.path)
        g_new = loadnpy('g_new')

        if PAR.SCHEME == 'GradientDescent':
            p_new = -g_new

        elif PAR.SCHEME == 'ConjugateGradient':
            # compute NLCG udpate
            p_new = cls.NLCG.compute()

        elif PAR.SCHEME =='QuasiNewton':
            # compute L-BFGS update
            if cls.iter == 1:
                p_new = -g_new
            else:
                cls.LBFGS.update()
                p_new = cls.LBFGS.solve()

                if cls.LBFGS.restarted:
                    cls.restart_search = 0

        else:
            raise ParameterError

        # save results
        unix.cd(cls.path)
        savenpy('p_new', p_new)
        savetxt('s_new', np.dot(g_new, p_new))


    ### line search methods

    # The following names are used in the line search methods:
    #     m - model
    #     p - search direction
    #     s - slope along search direction
    #     f - objective function value
    #     x - step length
    #     p_ratio - ratio of model norm to search direction norm
    #     s_ratio - ratio of current slope to previous slope

    def initialize_search(cls):
        """ Determines initial step length for line search
        """
        unix.cd(cls.path)

        m = loadnpy('m_new')
        p = loadnpy('p_new')
        f = loadtxt('f_new')
        norm_m = max(abs(m))
        norm_p = max(abs(p))
        p_ratio = float(norm_m/norm_p)

        if cls.iter > 1:
            s_ratio = loadtxt('s_new')/loadtxt('s_old')
            alpha = loadtxt('alpha')

        if not hasattr(cls, 'restart_search'):
            cls.restart_search = 0

        if cls.iter == 1:
            assert PAR.STEPLEN > 0.
            alpha = p_ratio * PAR.STEPLEN
        elif cls.restart_search:
            alpha *= 2.*s_ratio
        elif PAR.LINESEARCH in ['Bracket']:
            alpha *= 2.*s_ratio
        elif PAR.SCHEME in ['GradientDescent', 'ConjugateGradient']:
            alpha *= 2.*s_ratio
        else:
            alpha = 1.

        # ad hoc scaling
        if PAR.ADHOCSCALING:
            alpha *= PAR.ADHOCSCALING

        # limit maximum step length
        if PAR.STEPLENMAX:
            if alpha > p_ratio * PAR.STEPLENMAX:
                alpha = p_ratio * PAR.STEPLENMAX

        # reset search history
        cls.search_history = [[0., f]]
        cls.step_count = 0
        cls.restart_search = 0

        cls.isdone = 0
        cls.isbest = 0
        cls.isbrak = 0

        # write trial model
        savenpy('m_try', m + p*alpha)
        savetxt('alpha', alpha)

        cls.writer(cls.iter, 0., f)


    def search_status(cls):
        """ Determines status of line search

            Maintains line search history by keeping track of step length and
            function value from each trial model evaluation. From line search
            history, determines whether stopping criteria have been satisfied.
        """
        unix.cd(cls.path)

        # update search history
        x_ = loadtxt('alpha')
        f_ = loadtxt('f_try')

        if np.isnan(f_):
            raise ValueError

        cls.search_history += [[x_, f_]]
        cls.step_count += 1

        x = cls.step_lens()
        f = cls.func_vals()

        # is current step length the best so far?
        vals = cls.func_vals(sort=False)
        if np.all(vals[-1] < vals[:-1]):
            cls.isbest = 1

        # are stopping criteria satisfied?
        if PAR.LINESEARCH == 'Backtrack':
            if any(f[1:] < f[0]):
                cls.isdone = 1

        elif PAR.LINESEARCH == 'Bracket':
            if cls.isbrak:
                cls.isbest = 1
                cls.isdone = 1
            elif any(f[1:] < f[0]) and (f[-2] < f[-1]):
                cls.isbrak = 1

        elif PAR.LINESEARCH == 'Fixed':
            if any(f[1:] < f[0]) and (f[-2] < f[-1]):
                cls.isdone = 1

        # pass latest information to output writer
        cls.writer([], x_, f_)

        if cls.step_count >= PAR.SRCHMAX:
            cls.isdone = -1
            print ' line search failed [max iter]\n'

        return cls.isdone, cls.isbest


    def compute_step(cls):
        """ Computes next trial step length
        """
        unix.cd(cls.path)
        m = loadnpy('m_new')
        p = loadnpy('p_new')
        s = loadtxt('s_new')

        norm_m = max(abs(m))
        norm_p = max(abs(p))
        p_ratio = float(norm_m/norm_p)

        x = cls.step_lens()
        f = cls.func_vals()

        # compute trial step length
        if PAR.LINESEARCH == 'Backtrack':
            alpha = lib.backtrack2(f[0], s, x[1], f[1], b1=0.1, b2=0.5)

        elif PAR.LINESEARCH == 'Bracket':
            FACTOR = 2.
            if any(f[1:] < f[0]) and (f[-2] < f[-1]):
                alpha = lib.polyfit2(x, f)
            elif any(f[1:] < f[0]):
                alpha = loadtxt('alpha')*FACTOR
            else:
                alpha = loadtxt('alpha')*FACTOR**-1

        elif PAR.LINESEARCH == 'Fixed':
            alpha = p_ratio*(step + 1)*PAR.STEPLEN

        else:
            raise ValueError

        # write trial model
        savetxt('alpha', alpha)
        savenpy('m_try', m + p*alpha)


    def finalize_search(cls):
        """ Cleans working directory and writes updated model
        """
        unix.cd(cls.path)
        m = loadnpy('m_new')
        p = loadnpy('p_new')

        x = cls.step_lens()
        f = cls.func_vals()

        # clean working directory
        unix.rm('alpha')
        unix.rm('m_try')
        unix.rm('f_try')

        if cls.iter > 1:
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

        cls.writer([], [], [])


    ### line search utilities

    def step_lens(cls, sort=True):
        x, f = zip(*cls.search_history)
        x = np.array(x)
        f = np.array(f)
        f_sorted = f[abs(x).argsort()]
        x_sorted = x[abs(x).argsort()]
        if sort:
            return x_sorted
        else:
            return x


    def func_vals(cls, sort=True):
        x, f = zip(*cls.search_history)
        x = np.array(x)
        f = np.array(f)
        f_sorted = f[abs(x).argsort()]
        x_sorted = x[abs(x).argsort()]
        if sort:
            return f_sorted
        else:
            return f
