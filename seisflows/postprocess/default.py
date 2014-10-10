
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import exists
from seisflows.tools.configtools import getclass, ParameterObject

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

        self.path = PATH.POSTPROCESS
        unix.mkdir(self.path)


    def process_kernels(self,input=None,suffix='new'):
        """ Computes gradient and performs smoothing, preconditioning, and then
            scaling operations
        """
        assert(exists(input))
        unix.cd(input+'/'+'kernels')

        # combine kernels
        system.run( solver.combine,
            hosts='head',
            path=input+'/'+'kernels' )

        # write gradient
        g = solver.merge(solver.load('sum',type='kernel'))
        m = solver.merge(solver.load('../model',type='model'))
        g *= m
        solver.save(self.path+'/'+'gradient',solver.split(g))

        # apply smoothing
        if PAR.SMOOTH > 0.:
            system.run( solver.smooth,
                hosts='head',
                path=self.path,
                span=PAR.SMOOTH )

        # apply preconditioner
        if PATH.PRECOND:
            unix.cd(self.path)
            v = solver.merge(solver.load('gradient'))
            p = solver.merge(solver.load(PATH.PRECOND))
            unix.mv('gradient','gradient_noscale')
            solver.save('gradient',solver.split(v/p))

        # apply scaling
        if PAR.SCALE:
            unix.cd(self.path)
            g = solver.merge(solver.load('gradient',type='model'))
            g *= PAR.SCALE
            solver.save('gradient',solver.split(g))

