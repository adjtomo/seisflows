#
# This is Seisflows
#
# See LICENCE file
#
###############################################################################

# Import system modules
import sys
from os.path import join

# Import Numpy
import numpy as np

# Local imports
from seisflows.config import ParameterError
from seisflows.plugins import line_search, preconds
from seisflows.tools import msg, unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.math import angle
from seisflows.tools.seismic import Writer

# seisflows.config objects
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class base(object):
    """ Nonlinear optimization abstract base class

     Base class on top of which steepest descent, nonlinear conjugate, quasi-
     Newton and Newton methods can be implemented.  Includes methods for
     both search direction and line search.

     To reduce memory overhead, vectors are read from disk rather than passed
     from calling routines. For example, at the beginning of compute_direction
     the current gradient is  read from  'g_new' and the resulting search
     direction is written to 'p_new'. As the inversion progresses, other
     information is stored as well.

     Variables
        m_new - current model
        m_old - previous model
        m_try - line search model
        f_new - current objective function value
        f_old - previous objective function value
        f_try - line search function value
        g_new - current gradient direction
        g_old - previous gradient direction
        p_new - current search direction
        p_old - previous search direction
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # The default numerical parameters defined below should work well for a
        # range of applications without manual tuning. If the nonlinear
        # optimization procedure stagnates, it may be due to issues involving
        # data quality or the choice of data misfit, data processing, or
        # regularization parameters.  Problems in any of these areas usually
        # manifest themselves through stagnation of the nonlinear optimization
        # algorithm.

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
            step_count_max=PAR.STEPCOUNTMAX,
            path=PATH.WORKDIR+'/'+'output.optim')

        # prepare preconditioner
        if PAR.PRECOND:
            self.precond = getattr(preconds, PAR.PRECOND)()
        else:
            self.precond = None

        # prepare output logs
        self.writer = Writer(
                path=PATH.WORKDIR+'/'+'output.stats')

        # prepare scratch directory
        unix.mkdir(PATH.OPTIMIZE)
        if 'MODEL_INIT' in PATH:
            solver = sys.modules['seisflows_solver']
            self.save('m_new', solver.merge(solver.load(PATH.MODEL_INIT)))

    def compute_direction(self):
        """ Computes search direction
        """
        # the following code implements steepest descent
        # (for other algorithms, simply overload this method)
        g_new = self.load('g_new')
        if self.precond:
            p_new = -self.precond(g_new)
        else:
            p_new = -g_new
        self.save('p_new', p_new)

    def initialize_search(self):
        """ Determines first step length in line search
        """
        m = self.load('m_new')
        g = self.load('g_new')
        p = self.load('p_new')
        f = self.loadtxt('f_new')
        norm_m = max(abs(m))
        norm_p = max(abs(p))
        gtg = self.dot(g, g)
        gtp = self.dot(g, p)

        if self.restarted:
            self.line_search.clear_history()

        # optional step length safeguard
        if PAR.STEPLENMAX:
            self.line_search.step_len_max = \
                PAR.STEPLENMAX*norm_m/norm_p

        # determine initial step length
        alpha, _ = self.line_search.initialize(0., f, gtg, gtp)

        # optional initial step length override
        if PAR.STEPLENINIT and len(self.line_search.step_lens) <= 1:
            alpha = PAR.STEPLENINIT*norm_m/norm_p

        # write model corresponding to chosen step length
        self.savetxt('alpha', alpha)
        self.save('m_try', m + alpha*p)

    def update_search(self):
        """ Updates line search status and step length

          Status codes
              status > 0  : finished
              status == 0 : not finished
              status < 0  : failed
        """
        alpha, status = self.line_search.update(
            self.loadtxt('alpha'),
            self.loadtxt('f_try'))

        if status >= 0:
            # write model corresponding to chosen step length
            m = self.load('m_new')
            p = self.load('p_new')
            self.savetxt('alpha', alpha)
            self.save('m_try', m + alpha*p)
        return status

    def finalize_search(self):
        """ Prepares algorithm machinery and scratch directory for next
          model upate
        """
        m = self.load('m_new')
        g = self.load('g_new')
        p = self.load('p_new')
        x = self.line_search.search_history()[0]
        f = self.line_search.search_history()[1]

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
        self.writer('factor', -self.dot(g, g)**-0.5 * (f[1]-f[0])/(x[1]-x[0]))
        self.writer('gradient_norm_L1', np.linalg.norm(g, 1))
        self.writer('gradient_norm_L2', np.linalg.norm(g, 2))
        self.writer('misfit', f[0])
        self.writer('restarted', self.restarted)
        self.writer('slope', (f[1]-f[0])/(x[1]-x[0]))
        self.writer('step_count', self.line_search.step_count)
        self.writer('step_length', x[f.argmin()])
        self.writer('theta', 180.*np.pi**-1*angle(p, -g))

        self.line_search.writer.newline()

    def retry_status(self):
        """ Determines if restart is worthwhile

          After failed line search, determines if restart is worthwhile by
          checking, in effect, if search direction was the same as gradient
          direction
        """
        g = self.load('g_new')
        p = self.load('p_new')
        theta = angle(p, -g)

        if PAR.VERBOSE >= 2:
            print ' theta: %6.3f' % theta

        thresh = 1.e-3
        if abs(theta) < thresh:
            return 0
        else:
            return 1

    def restart(self):
        """ Restarts nonlinear optimization algorithm

          Keeps current position in model space, but discards history of
          nonlinear optimization algorithm in an attempt to recover from
          numerical stagnation
        """
        g = self.load('g_new')
        self.save('p_new', -g)
        self.line_search.clear_history()
        self.restarted = 1
        self.line_search.writer.iter -= 1
        self.line_search.writer.newline()

    def dot(self, x, y):
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
        # writes scalars to disk
        np.savetxt(PATH.OPTIMIZE+'/'+filename, [scalar], '%11.6e')
