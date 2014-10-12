#!/usr/bin/env python

import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import abspath, join, loadtxt, savetxt
from seisflows.tools.configtools import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

from problems import rosenbrock as problem



class run(getclass('workflow','inversion')):
  """ Optimization unit test.

      Tests nonlinear optimization procedure using inexpensive test function.
  """

  def __init__(cls):

    # check paths
    if 'GLOBAL' not in PATH:
       setattr(PATH,'GLOBAL',abspath('./scratch'))

    if 'SUBMIT' not in PATH:
       setattr(PATH,'SUBMIT',abspath('.'))

    if 'OPTIMIZE' not in PATH:
       setattr(PATH,'OPTIMIZE',PATH.GLOBAL)


    # check parameters
    if 'BEGIN' not in PAR:
        raise Exception

    if 'END' not in PAR:
        raise Exception


  def setup(cls):
    cls.path = PATH.GLOBAL

    unix.mkdir(cls.path)
    unix.cd(cls.path)

    cls.optimize = getclass('optimize','default')()
    cls.optimize.setup()

    # prepare starting model
    unix.cd(cls.path)
    m = np.array([-1.2,1])
    savenpy('m_new',m)


  def initialize(cls):
    pass 


  def search_status(cls):
    cls.evaluate_function()
    isdone, _ = cls.optimize.search_status()
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
    m = loadnpy('m_lcg')
    if PAR.SCHEME == 'GaussNewton':
      pass
    elif PAR.SCHEME == 'Newton':
      g = problem.grad(m)
      savenpy('g_lcg',g)


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

