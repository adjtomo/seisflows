
from os.path import join
import sys
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists, loadtxt, savetxt
from seisflows.config import ParameterError 

from seisflows.tools.math import angle, polyfit2, backtrack2
from seisflows.optimize.lib.io import Writer, StepWriter


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class steepest_descent(object):
    """ Steepest descent

     Default numerical parameters provided below should work well for a wide range
     inversions without the need for manual tuning.

     To reduce memory overhead, vectors are read from disk rather than passed
     from a calling routine. At the start of each search direction computation
     the current model and gradient are read from files 'm_new' and 'g_new';
     the resulting search direction is written to 'p_new'. As the inversion
     progresses, other information is stored to disk as well.
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        if 'SUBMIT' not in PATH:
            raise ParameterError

        if 'OPTIMIZE' not in PATH:
            setattr(PATH, 'OPTIMIZE', join(PATH.SCRATCH, 'optimize'))

        if 'MODEL_INIT' not in PATH:
            setattr(PATH, 'MODEL_INIT', None)

        # preconditioner
        if 'PRECOND' not in PAR:
            setattr(PAR, 'PRECOND', None)

        # line search algorithm
        if 'LINESEARCH' not in PAR:
            if  PAR.OPTIMIZE in ['LBFGS']:
                setattr(PAR, 'LINESEARCH', 'Backtrack')
            else:
                setattr(PAR, 'LINESEARCH', 'Bracket')

        # maximum number of trial steps
        if 'STEPMAX' not in PAR:
            setattr(PAR, 'STEPMAX', 10)

        # initial step length as percent of current model
        if 'STEPINIT' not in PAR:
            setattr(PAR, 'STEPINIT', 0.05)

        # optional initial step length safegaurd
        if 'STEPTHRESH' not in PAR:
            setattr(PAR, 'STEPTHRESH', None)

        # step length factor in bracketing line search
        if 'STEPFACTOR' not in PAR:
            setattr(PAR, 'STEPFACTOR', 0.5)

        # optional parameter, can be useful for NLCG line search
        if 'STEPOVERSHOOT' not in PAR:
            setattr(PAR, 'STEPOVERSHOOT', 0.)

        # ad hoc factor by which to scale gradient
        if 'ADHOCFACTOR' not in PAR:
            setattr(PAR, 'ADHOCFACTOR', 1.)


    def setup(self):
        """ Sets up nonlinear optimization machinery
        """
        # prepare output writers
        self.writer = Writer(
                path=PATH.OUTPUT)

        self.stepwriter = StepWriter(
                path=PATH.SUBMIT)

        # write initial model
        if exists(PATH.MODEL_INIT):
            solver = sys.modules['seisflows_solver']
            src = PATH.MODEL_INIT
            dst = join(PATH.OPTIMIZE, 'm_new')
            savenpy(dst, solver.merge(solver.load(src)))


    def precond(self):
        """ Loads preconditioner machinery
        """
        from seisflows.plugins import preconds

        if PAR.PRECOND in dir(preconds):
            return getattr(preconds, PAR.PRECOND)()
        elif PAR.PRECOND:
            return getattr(preconds, 'diagonal')()
        else:
            return None


    # The following names are used in the 'compute_direction' method and for
    # writing information to disk:
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
        g_new = self.load('g_new')
        p_new, self.restarted = -g_new, False
        self.save('p_new', p_new)
        savetxt('s_new', self.dot(g_new, p_new))
        return p_new


    # The following names are used exclusively for the line search:
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
        m = self.load('m_new')
        p = self.load('p_new')
        f = self.loadtxt('f_new')
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
            alpha = p_ratio*PAR.STEPINIT
        elif self.restarted:
            alpha = p_ratio*PAR.STEPINIT
        elif PAR.OPTIMIZE in ['LBFGS']:
            alpha = 1.
        else:
            alpha = self.initial_step()

        # optional ad hoc scaling
        if PAR.STEPOVERSHOOT:
            alpha *= PAR.STEPOVERSHOOT

        # optional maximum step length safegaurd
        if PAR.STEPTHRESH:
            if alpha > p_ratio * PAR.STEPTHRESH and \
                self.iter > 1:
                alpha = p_ratio * PAR.STEPTHRESH

        # write trial model corresponding to chosen step length
        savetxt('alpha', alpha)
        self.save('m_try', m + alpha*p)

        # append latest statistics
        self.stepwriter(steplen=0., funcval=f)


    def update_status(self):
        """ Updates line search status

            Maintains line search history by keeping track of step length and
            function value from each trial model evaluation. From line search
            history, determines whether stopping criteria have been satisfied.
        """
        x_ = self.loadtxt('alpha')
        f_ = self.loadtxt('f_try')
        if np.isnan(f_):
            raise ValueError

        # update search history
        self.search_history += [[x_, f_]]
        self.step_count += 1
        x = self.step_lens()
        f = self.func_vals()

        fmin = f.min()
        imin = f.argmin()

        # is current step length the best so far?
        vals = self.func_vals(sort=False)
        if np.all(vals[-1] < vals[:-1]):
            self.isbest = 1

        # are stopping criteria satisfied?
        if PAR.LINESEARCH == 'Fixed':
            if (fmin < f[0]) and any(fmin < f[imin:]):
                self.isdone = 1

        elif PAR.LINESEARCH == 'Bracket' or \
            self.iter == 1 or self.restarted:
            if self.isbrak:
                self.isdone = 1
            elif (fmin < f[0]) and any(fmin < f[imin:]):
                self.isbrak = 1

        elif PAR.LINESEARCH == 'Backtrack':
            if fmin < f[0]:
                self.isdone = 1

        # append latest statistics
        self.stepwriter(steplen=x_, funcval=f_)

        return self.isdone


    def compute_step(self):
        """ Computes next trial step length
        """
        m = self.load('m_new')
        g = self.load('g_new')
        p = self.load('p_new')
        s = self.loadtxt('s_new')

        norm_m = max(abs(m))
        norm_p = max(abs(p))
        p_ratio = float(norm_m/norm_p)

        x = self.step_lens()
        f = self.func_vals()

        # compute trial step length
        if PAR.LINESEARCH == 'Fixed':
            alpha = p_ratio*(self.step_count + 1)*PAR.STEPINIT

        elif PAR.LINESEARCH == 'Bracket' or \
            self.iter==1 or self.restarted:
            if any(f[1:] < f[0]) and (f[-2] < f[-1]):
                alpha = polyfit2(x, f)

            elif any(f[1:] <= f[0]):
                alpha = self.loadtxt('alpha')*PAR.STEPFACTOR**-1
            else:
                alpha = self.loadtxt('alpha')*PAR.STEPFACTOR

        elif PAR.LINESEARCH == 'Backtrack':
            # calculate slope along 1D profile
            slope = s/self.dot(g,g)**0.5
            if PAR.ADHOCFACTOR:
                slope *= PAR.ADHOCFACTOR            

            alpha = backtrack2(f[0], slope, x[1], f[1], b1=0.1, b2=0.5)

        # write trial model corresponding to chosen step length
        self.savetxt('alpha', alpha)
        self.save('m_try', m + alpha*p)


    def finalize_search(self):
        """ Cleans working directory and writes updated model
        """
        m = self.load('m_new')
        g = self.load('g_new')
        p = self.load('p_new')
        s = self.loadtxt('s_new')

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
        self.save('m_new', m + alpha*p)
        savetxt('f_new', f.min())

        # append latest statistics
        self.writer('factor', -self.dot(g,g)**-0.5 * (f[1]-f[0])/(x[1]-x[0]))
        self.writer('gradient_norm_L1', np.linalg.norm(g, 1))
        self.writer('gradient_norm_L2', np.linalg.norm(g, 2))
        self.writer('misfit', f[0])
        self.writer('restarted', self.restarted)
        self.writer('slope', (f[1]-f[0])/(x[1]-x[0]))
        self.writer('step_count', self.step_count)
        self.writer('step_length', x[f.argmin()])
        self.writer('theta', 180.*np.pi**-1*angle(p,-g))

        self.stepwriter.newline()


    def retry_status(self):
        """ Returns false if search direction was the same as gradient
          direction; returns true otherwise
        """
        g = self.load('g_new')
        p = self.load('p_new')

        thresh = 1.e-3
        theta = angle(p,-g)

        if PAR.VERBOSE >= 2:
            print ' theta: %6.3f' % theta

        if abs(theta) < thresh:
            return 0
        else:
            return 1


    def restart(self):
        """ Discards history of algorithm; prepares to start again from 
          gradient direction
        """
        g = self.load('g_new')
        self.save('p_new', -g)
        savetxt('s_new', self.dot(g,g))
        self.restarted = 1
        self.stepwriter.iter -= 1
        self.stepwriter.newline()



    ### line search utilities

    def initial_step(self):
        """ Determines first trial step in line search; see eg Nocedal and 
          Wright 2e section 3.5
        """
        alpha = self.loadtxt('alpha')
        s_new = self.loadtxt('s_new')
        s_old = self.loadtxt('s_old')
        s_ratio = s_new/s_old
        return 2.*s_ratio*alpha


    def step_lens(self, sort=True):
        """ Returns previous step lengths from search history
        """
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
        """ Returns previous function values from search history
        """
        x, f = zip(*self.search_history)
        x = np.array(x)
        f = np.array(f)
        f_sorted = f[abs(x).argsort()]
        x_sorted = x[abs(x).argsort()]
        if sort:
            return f_sorted
        else:
            return f


    ### utilities

    def dot(self,x,y):
        """ Computes inner product between vectors
        """
        return np.dot(
            np.squeeze(x),
            np.squeeze(y))

    def load(self, filename):
        return loadnpy(PATH.OPTIMIZE+'/'+filename)

    def save(self, filename, v):
        savenpy(PATH.OPTIMIZE+'/'+filename, v)


    def loadtxt(self, filename):
        return loadtxt(PATH.OPTIMIZE+'/'+filename)

    def savetxt(self, filename, c):
        savetxt(PATH.OPTIMIZE+'/'+filename, c)


