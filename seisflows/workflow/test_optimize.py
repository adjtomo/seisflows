
import sys
import numpy as np

from os.path import abspath
from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import savetxt
from seisflows.config import ParameterError

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

optimize = sys.modules['seisflows_optimize']

import problems.rosenbrock as problem


class test_optimize(object):
    """ Optimization unit test.

        Tests nonlinear optimization procedure with inexpensive test function.
    """

    def check(cls):
        cls.path = PATH.OPTIMIZE

        # check parameters
        if 'OPTIMIZE' not in PAR: 
            setattr(PAR, 'OPTIMIZE', 'base')

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        if 'BEGIN' not in PAR:
            raise Exception

        if 'END' not in PAR:
            raise Exception

        if 'THRESH' not in PAR:
            setattr(PAR, 'THRESH', 1.e-5)

        # check paths
        if 'SCRATCH' not in PATH:
            setattr(PATH,'SCRATCH',abspath('./scratch'))

        if 'SUBMIT' not in PATH:
            setattr(PATH,'SUBMIT',abspath('.'))

        if 'OPTIMIZE' not in PATH:
            setattr(PATH,'OPTIMIZE',PATH.SCRATCH)

        # assertions
        assert PAR.BEGIN == 1
        assert PAR.BEGIN < PAR.END

        assert 'func' in dir(problem)
        assert 'grad' in dir(problem)
        assert 'model_init' in dir(problem)


    def main(cls):
        cls.setup()

        for cls.iter in range(PAR.BEGIN, PAR.END+1):
            print 'Starting iteration', cls.iter
            optimize.iter = cls.iter

            print "Computing search direction"
            cls.compute_direction()

            print "Computing step length"
            cls.line_search()

            cls.finalize()
            print ''


    def setup(cls):
        unix.mkdir(cls.path)
        unix.cd(cls.path)

        optimize.check()
        optimize.setup()

        unix.cd(cls.path)
        m = problem.model_init()
        savenpy('m_new',m)


    def compute_direction(cls):
        cls.evaluate_gradient()
        if PAR.OPTIMIZE in ['Newton', 'GaussNewton']:
            cls.compute_direction_newton()
        else:
            optimize.compute_direction()


    def compute_direction_newton(cls):
        optimize.initialize_newton()

        for ilcg in range(PAR.LCGMAX):
            m = loadnpy('m_lcg')
            g = problem.grad(m)
            savenpy('g_lcg', g)
            isdone = optimize.iterate_newton()
            if isdone:
                break


    def line_search(cls):
        optimize.initialize_search()

        while True:
            cls.evaluate_function()
            optimize.update_status()

            if optimize.isdone:
                optimize.finalize_search()
                break

            elif optimize.step_count < PAR.STEPMAX:
                optimize.compute_step()
                continue

            else:
                retry = optimize.retry_status()
                if retry:
                    print ' Line search failed... retry'
                    optimize.restart()
                    cls.line_search()
                    break
                else:
                    print ' Line search failed... abort'
                    sys.exit(-1)


    def evaluate_function(cls):
        m = loadnpy('m_try')
        f = problem.func(m)
        savetxt('f_try',f)


    def evaluate_gradient(cls):
        m = loadnpy('m_new')
        f = problem.func(m)
        g = problem.grad(m)
        savetxt('f_new',f)
        savenpy('g_new',g)


    def finalize(cls):
        m_new = loadnpy('m_new')
        m_old = loadnpy('m_old')

        if PAR.VERBOSE > 0:
            print '%14.7e %14.7e'%tuple(m_new)

        if cls.status(m_new, m_old):
            print 'Stopping criteria met.\n'
            np.savetxt('niter', [cls.iter], '%d')
            sys.exit(0)

        elif cls.iter >= PAR.END:
            print 'Maximum number of iterations exceeded.\n'
            sys.exit(-1)


    def status(cls, m_new, m_old):
        _finished = True
        if np.linalg.norm(m_new-m_old)/np.linalg.norm(m_new) < PAR.THRESH:
            return _finished
        else:
            return not _finished

