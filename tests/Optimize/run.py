#!/usr/bin/env python

import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import abspath, irange, join, loadtxt, savetxt
from seisflows.tools.configtools import loadclass, loadvars, ParameterObj

from problems import rosenbrock as problem

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

PAR.update(loadvars('parameters','.'))
PATH.update(loadvars('paths','.'))


class run(object):
    """ Optimization unit test.

        Tests nonlinear optimization procedure using inexpensive test function.
    """

    def __init__(cls):

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


        # instantiate
        global optimize
        optimize = loadclass('optimize',PAR.OPTIMIZE)()



    def main(cls):
        cls.setup()

        for cls.iter in irange(PAR.BEGIN,PAR.END):
            print 'Starting iteration', cls.iter
            optimize.iter = cls.iter

            print "Computing search direction"
            cls.compute_direction()

            print "Computing step length"
            cls.line_search()

            cls.finalize()
            print ''

            if cls.isdone:
                return


    def setup(cls):
        unix.mkdir(cls.path)
        unix.cd(cls.path)

        optimize.check()
        optimize.setup()

        # prepare starting model
        unix.cd(cls.path)
        m = np.array([-1.2,1])
        savenpy('m_new',m)


    def compute_direction(cls):
        cls.evaluate_gradient()
        optimize.compute_direction()


    def line_search(cls):
        optimize.initialize_search()

        for cls.step in irange(1,PAR.SRCHMAX):
            isdone = cls.search_status()

            if isdone==1:
                optimize.finalize_search()
                break
            elif isdone==0:
                optimize.compute_step()
                continue
            elif isdone==-1:
                cls.isdone = -1
                print ' line search failed'


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


    def apply_hessian(cls):
        raise NotImplementedError


    def finalize(cls):
        print '%14.7e %14.7e'%tuple(loadnpy('m_new'))

        m_new = loadnpy('m_new')
        m_old = loadnpy('m_old')

        # check stopping condition
        d = np.linalg.norm(m_new-m_old)/np.linalg.norm(m_new)
        if d < 1.e-6:
            cls.isdone = True
            cls.msg = 'Stopping criteria met.\n'
        else:
            cls.isdone = False


# run
if __name__ == '__main__':
    run().main()
