
import sys
import numpy as np

from os.path import join
from seisflows.config import ParameterError
from seisflows.plugins import line_search, preconds
from seisflows.tools import msg, unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.math import angle, polyfit2, backtrack2
from seisflows.tools.seismic import  Writer, StepWriter


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class base(object):
    """ Abstract base class

     Default numerical parameters provided below should work well for a wide range
     inversions without the need for manual tuning. If the nonlinear optimization
     procedure stagnates, it may be due to the objective function rather than the 
     numerical parameters.

     To reduce memory overhead, vectors are read from disk rather than passed
     from a calling routine. At the start of each search direction computation
     the current model and gradient are read from files 'm_new' and 'g_new';
     the resulting search direction is written to 'p_new'. As the inversion
     progresses, other information is stored to disk as well.
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # line search algorithm
        if 'LINESEARCH' not in PAR:
            setattr(PAR, 'LINESEARCH', 'Bracket')

        # preconditioner
        if 'PRECOND' not in PAR:
            setattr(PAR, 'PRECOND', None)

        # maximum number of trial steps
        if 'STEPCOUNTMAX' not in PAR:
            setattr(PAR, 'STEPCOUNTMAX', 10)

        # initial step length as fraction of current model
        if 'STEPLENINIT' not in PAR:
            setattr(PAR, 'STEPLENINIT', 0.05)

        # maximum step length as fraction of current model
        if 'STEPLENMAX' not in PAR:
            setattr(PAR, 'STEPLENMAX', 0.5)

        # where temporary files are written
        if 'OPTIMIZE' not in PATH:
            setattr(PATH, 'OPTIMIZE', PATH.SCRATCH+'/'+'optimize')


        # assertions
        if 'WORKDIR' not in PATH:
            raise ParameterError

        if PAR.OPTIMIZE in ['base']:
            print msg.CompatibilityError1
            sys.exit(-1)

        if PAR.LINESEARCH:
            assert PAR.LINESEARCH in dir(line_search)

        if PAR.PRECOND:
            assert PAR.PRECOND in dir(preconds)

        if PAR.STEPLENINIT:
            assert 0. < PAR.STEPLENINIT

        if PAR.STEPLENMAX:
            assert 0. < PAR.STEPLENMAX

        if PAR.STEPLENINIT and PAR.STEPLENMAX:
            assert PAR.STEPLENINIT < PAR.STEPLENMAX


    def setup(self):
        """ Sets up nonlinear optimization machinery
        """
        # prepare line search machinery
        self.line_search = getattr(line_search, PAR.LINESEARCH)(
            step_count_max=PAR.STEPCOUNTMAX)

        if PAR.PRECOND:
            self.precond = getattr(preconds, PAR.PRECOND)()
        else:
            self.precond = None

        # prepare output writers
        self.writer = Writer(
                path=PATH.WORKDIR+'/'+'output.stats')

        self.stepwriter = StepWriter(
                path=PATH.WORKDIR+'/'+'output.optim')

        # prepare scratch directory
        unix.mkdir(PATH.OPTIMIZE)
        if 'MODEL_INIT' in PATH:
            solver = sys.modules['seisflows_solver']
            self.save('m_new', solver.merge(solver.load(PATH.MODEL_INIT)))


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

    def compute_direction(self):
        """ Computes model update direction from stored gradient
        """
        # must be implemented by subclass
        raise NotImplementedError


    def initialize_search(self):
        """ Determines first step length in line search
        """
        m = self.load('m_new')
        g = self.load('g_new')
        p = self.load('p_new')
        f = self.loadtxt('f_new')
        norm_m = max(abs(m))
        norm_p = max(abs(p))

        self.line_search.step_count = 0
        if self.restarted: self.line_search.clear_history()
        self.line_search.step_lens += [0.]
        self.line_search.func_vals += [f]
        self.line_search.gtg += [self.dot(g,g)]
        self.line_search.gtp += [self.dot(g,p)]

        self.stepwriter(
            steplen=0., 
            funcval=self.loadtxt('f_new'))

        # determine initial step length
        if PAR.STEPLENMAX:
            self.line_search.step_len_max = \
                PAR.STEPLENMAX*norm_m/norm_p

        if PAR.STEPLENINIT and len(self.line_search.step_lens)<=1:
            alpha = PAR.STEPLENINIT*norm_m/norm_p

        else:
            alpha = self.line_search.initial_step()

        # write model corresponding to chosen step length
        self.savetxt('alpha', alpha)
        self.save('m_try', m + alpha*p)


    def update_search(self):
        """ Updates line search status
        """
        self.line_search.step_count += 1
        self.line_search.step_lens += [self.loadtxt('alpha')]
        self.line_search.func_vals += [self.loadtxt('f_try')]

        self.stepwriter(
            steplen=self.loadtxt('alpha'),
            funcval=self.loadtxt('f_try'))

        alpha, status = self.line_search.update()
        if status >= 0:
            # write model corresponding to chosen step length
            m = self.load('m_new')
            p = self.load('p_new')
            self.savetxt('alpha', alpha)
            self.save('m_try', m + alpha*p)
        return status


    def finalize_search(self):
        """ Writes output statistics and prepares scratch directory for next
          model upate
        """
        m = self.load('m_new')
        g = self.load('g_new')
        p = self.load('p_new')
        x = self.line_search.current_vals()[0]
        f = self.line_search.current_vals()[1]

        # clean scratch directory
        unix.cd(PATH.OPTIMIZE)
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

        unix.mv('m_try', 'm_new')
        self.savetxt('f_new', f.min())

        # output latest statistics
        self.writer('factor', -self.dot(g,g)**-0.5 * (f[1]-f[0])/(x[1]-x[0]))
        self.writer('gradient_norm_L1', np.linalg.norm(g, 1))
        self.writer('gradient_norm_L2', np.linalg.norm(g, 2))
        self.writer('misfit', f[0])
        self.writer('restarted', self.restarted)
        self.writer('slope', (f[1]-f[0])/(x[1]-x[0]))
        self.writer('step_count', self.line_search.step_count)
        self.writer('step_length', x[f.argmin()])
        self.writer('theta', 180.*np.pi**-1*angle(p,-g))
        self.stepwriter.newline()


    def retry_status(self):
        """ Determines if retry is worthwhile after failed line search

          Determines if retry is worthwhile by checking, in effect, if search 
          direction was the same as gradient direction
        """
        g = self.load('g_new')
        p = self.load('p_new')
        theta = angle(p,-g)

        if PAR.VERBOSE >= 2:
            print ' theta: %6.3f' % theta

        thresh = 1.e-3
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
        self.line_search.clear_history()
        self.line_search.step_count = 0
        self.restarted = 1
        self.stepwriter.iter -= 1
        self.stepwriter.newline()


    def dot(self,x,y):
        """ Computes inner product between vectors
        """
        return np.dot(
            np.squeeze(x),
            np.squeeze(y))

    def load(self, filename):
        # reads vectors from disk
        return loadnpy(PATH.OPTIMIZE+'/'+filename)

    def save(self, filename, array):
        # writes vectors to disk
        savenpy(PATH.OPTIMIZE+'/'+filename, array)

    def loadtxt(self, filename):
        # reads scalars from disk
        return float(np.loadtxt(PATH.OPTIMIZE+'/'+filename))

    def savetxt(self, filename, scalar):
        # writes scalars
        np.savetxt(PATH.OPTIMIZE+'/'+filename, [scalar], '%11.6e')


