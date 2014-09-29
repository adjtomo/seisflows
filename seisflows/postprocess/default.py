
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.configure import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()


class default(object):
  """ Postprocessing class

    First, combines kernels (i.e. contributions from individual sources) to 
    obtain the gradient direction. Next, performs smoothing, preconditioning,
    and scaling operations on gradient in accordance with parameter settings.
  """

  def __init__(self):
    """ Constructor
    """
    # check user supplied parameters
    if 'SMOOTH' not in PAR:
	setattr(PAR,'SMOOTH',0.)

    if 'SCALE' not in PAR:
	setattr(PAR,'SCALE',1.)


  def process_kernels(self):
    # combine kernels
    system.run( solver.combine, 
	hosts='head',
	path=PATH.GRAD+'/'+'kernels' )

    # construct mask
    unix.cd(PATH.GRAD+'/'+'kernels')
    g = solver.merge(solver.load('sum',type='kernel'))
    m = solver.merge(solver.load('../model',type='model'))
    mask = m>0

    # write gradient
    g[mask] = g[mask]/m[mask]
    g[np.invert(mask)] = 0.
    solver.save(PATH.GRAD+'/'+'grad',solver.split(g))

    # apply smoothing
    if PAR.SMOOTH > 0.:
      system.run( solver.smooth, 
	  hosts='head',
	  path=PATH.GRAD,
	  span=PAR.SMOOTH )

    # apply preconditioner
    if PATH.PRECOND:
      unix.cd(PATH.GRAD)
      v = solver.merge(solver.load('grad'))
      p = solver.merge(solver.load(PATH.PRECOND))
      unix.mv('grad','grad_noscale')
      solver.save('grad',solver.split(v/p))

    # apply scaling
    if PAR.SCALE:
      unix.cd(PATH.GRAD)
      g = solver.merge(solver.load('grad',type='model'))
      g *= PAR.SCALE

    savenpy(PATH.OPTIMIZE+'/'+'g_new',g)


