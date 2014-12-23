import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import solver


class ChenTromp(object):
    """ Postprocessing class

      Combines contributions from individual sources to obtain the gradient
      direction, and performs scaling, clipping, smoothing, and preconditioning
      operations on gradient in accordance with parameter settings.
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # check postprocessing settings
        if 'SCALE' not in PAR:
            setattr(PAR, 'SCALE', False)

        if 'CLIP' not in PAR:
            setattr(PAR, 'CLIP', 0.)

        if 'SMOOTH' not in PAR:
            setattr(PAR, 'SMOOTH', 0.)

        if 'PRECOND' not in PATH:
            setattr(PATH, 'PRECOND', None)


    def setup(self):
        # nothing to do
        pass


    def merge_kernels(self, tag, path):
        """ Converts kernel dictionary to gradient vector
        """
        unix.cd(path + '/' + 'kernels')
        g = solver.merge(solver.load('sum', type='kernel', verbose=True))
        return g


    def process_kernels(self, tag='grad', path=None):
        """ Computes gradient and performs scaling, smoothing, and 
          preconditioning operations
        """
        assert (exists(path))

        # combine kernels
        system.run('solver', 'combine',
                   hosts='head',
                   path=path + '/' + 'kernels')

        g = self.merge_kernels(tag, path)

        # apply ad hoc scaling
        if float(PAR.SCALE) != 1.: g *= PAR.SCALE
        solver.save(path + '/' + tag, solver.split(g))

        # apply clipping
        if PAR.CLIP > 0.:
            raise NotImplementedError

        # apply smoothing
        if PAR.SMOOTH > 0.:
            system.run('solver', 'smooth',
                       hosts='head',
                       path=path,
                       span=PAR.SMOOTH)

        # apply preconditioner
        if PATH.PRECOND:
            unix.cd(path)
            g /= solver.merge(solver.load(PATH.PRECOND))
            unix.mv(tag, tag + '_noprecond')
            solver.save(tag, solver.split(g))

        if 'OPTIMIZE' in PATH:
            if tag == 'grad':
                g = solver.merge(solver.load(path +'/'+ tag, type='model', verbose=True))
                savenpy(PATH.OPTIMIZE +'/'+ 'g_new', g)

