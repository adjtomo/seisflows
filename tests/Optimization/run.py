#!/usr/bin/python

import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import join, loadtxt, savetxt
from seisflows.tools.configure import getclass, ParameterObject

import parameters
import paths

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')


from problems import rosenbrock as problem



class run(getclass('workflow','inversion')):
  """ Optimization unit run.

      Tests nonlinear optimization procedure with inexpensive test function.
  """

  def __init__(self):
    pass


  def setup(self):
    cwd = join(unix.pwd(),'scratch')
    unix.mkdir(cwd)
    unix.cd(cwd)

    self.optimize = getclass('optimize','default')(
      path=cwd,
      output=cwd+'/'+'output.optim' )

    # prepare starting model
    unix.cd(cwd)
    m = np.array([-1.2,1])
    savenpy('m_new',m)


  def initialize(self):
    pass 


  def search_status(self):
    self.evaluate_function()
    isdone, _ = self.optimize.search_status()
    return isdone


  def evaluate_function(self):
    m = loadnpy('m_try')
    f = problem.func(m)
    savetxt('f_try',f)


  def evaluate_gradient(self):
    m = loadnpy('m_new')
    f = problem.func(m)
    g = problem.grad(m)
    savetxt('f_new',f)
    savenpy('g_new',g)


  def apply_hessian(self):
    m = loadnpy('m_lcg')
    if PAR.SCHEME == 'gn':
      pass
    elif PAR.SCHEME == 'tn':
      g = problem.grad(m)
      savenpy('g_lcg',g)


  def finalize(self):
    print '%14.7e %14.7e'%tuple(loadnpy('m_new'))

    m_new = loadnpy('m_new')
    m_old = loadnpy('m_old')

    # check stopping condition
    d = np.linalg.norm(m_new-m_old)/np.linalg.norm(m_new)
    if d < 1.e-6:
      self.isdone = True
      self.msg = 'Stopping criteria met.\n'
    else:
      self.isdone = False


# run run
if __name__ == '__main__':
    run().main()

