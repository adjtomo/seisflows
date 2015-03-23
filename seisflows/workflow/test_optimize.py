
from os.path import abspath
import sys

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import savetxt
from seisflows.tools.config import ParameterObj
import problems.rosenbrock as problem

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import optimize


class test_optimize(object):
    """ Optimization unit test.

        Tests nonlinear optimization procedure with inexpensive test function.
    """

    def check(cls):

        # check parameters
        if 'OPTIMIZE' not in PAR: 
            setattr(PAR,'OPTIMIZE','default')

        if 'BEGIN' not in PAR:
            raise Exception

        if 'END' not in PAR:
            raise Exception


        # check paths
        if 'GLOBAL' not in PATH:
            setattr(PATH,'GLOBAL',abspath('./scratch'))

        if 'SUBMIT' not in PATH:
            setattr(PATH,'SUBMIT',abspath('.'))

        if 'OPTIMIZE' not in PATH:
            setattr(PATH,'OPTIMIZE',PATH.GLOBAL)

        cls.path = PATH.GLOBAL


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
        m = np.array([-1.2,1])
        savenpy('m_new',m)


    def compute_direction(cls):
        cls.evaluate_gradient()
        if PAR.OPTIMIZE == 'Newton':
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

        for cls.step in range(1, PAR.SRCHMAX+1):
            isdone = cls.search_status()
            if isdone==1:
                optimize.finalize_search()
                break
            elif isdone==0:
                optimize.compute_step()
                continue
            elif isdone==-1:
                sys.exit()


    def search_status(cls):
        cls.evaluate_function()
        isdone, _ = optimize.search_status()
        return isdone


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
        print '%14.7e %14.7e'%tuple(loadnpy('m_new'))

        m_new = loadnpy('m_new')
        m_old = loadnpy('m_old')

        d = np.linalg.norm(m_new-m_old)/np.linalg.norm(m_new)
        if d < 1.e-5:
            print 'Stopping criteria met.\n'
            sys.exit()

