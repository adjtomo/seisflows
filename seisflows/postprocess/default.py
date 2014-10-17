
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import exists
from seisflows.tools.configtools import getclass, GlobalStruct

PAR = GlobalStruct('parameters')
PATH = GlobalStruct('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()


class default(object):
    """ Postprocessing class

      First, combines kernels (i.e. contributions from individual sources) to
      obtain the gradient direction. Next, performs smoothing, preconditioning,
      and scaling operations on gradient in accordance with parameter settings.
    """

    def __init__(self):

        # check postprocessing settings
        if 'PRECOND' not in PATH:
            setattr(PATH,'PRECOND',None)

        if 'SMOOTH' not in PAR:
            setattr(PAR,'SMOOTH',0.)

        if 'SCALE' not in PAR:
            setattr(PAR,'SCALE',1.)


    def process_kernels(self,tag='grad',path=None,optim_path=None):
        """ Computes gradient and performs smoothing, preconditioning, and then
            scaling operations
        """
        assert(exists(path))

        # combine kernels
        system.run( solver.combine,
            hosts='head',
            path=path+'/'+'kernels' )

        # write gradient
        unix.cd(path+'/'+'kernels')
        g = solver.merge(solver.load('sum',type='kernel'))
        g *= solver.merge(solver.load('../model',type='model'))
        solver.save(path+'/'+tag,solver.split(g))

        # apply smoothing
        if PAR.SMOOTH > 0.:
            system.run( solver.smooth,
                hosts='head',
                path=path,
                span=PAR.SMOOTH )

        # apply preconditioner
        if PATH.PRECOND:
            unix.cd(path)
            g = solver.merge(solver.load(tag))
            p = solver.merge(solver.load(PATH.PRECOND))
            unix.mv(tag,tag+'_noscale')
            solver.save(tag,solver.split(g/p))

        # apply scaling
        if PAR.SCALE:
            unix.cd(path)
            g = solver.merge(solver.load(tag,type='model'))
            g *= PAR.SCALE
            solver.save(tag,solver.split(g))

        if optim_path:
            savenpy(optim_path,g)

