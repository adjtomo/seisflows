
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class base(object):
    """ Postprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        # check postprocessing settings
        if 'CLIP' not in PAR:
            setattr(PAR, 'CLIP', 0.)

        if 'SMOOTH' not in PAR:
            setattr(PAR, 'SMOOTH', 0.)

        if 'REGULARIZE' not in PAR:
            setattr(PAR, 'REGULARIZE', False)

        if 'PRECOND' not in PATH:
            setattr(PATH, 'PRECOND', None)

        if 'FACTOR' not in PAR:
            setattr(PAR, 'FACTOR', False)

        if PATH.PRECOND:
            assert exists(PATH.PRECOND)


    def setup(self):
        """ Performs any required initialization or setup tasks
        """
        pass


    def write_gradient(self, path):
        """ Writes gradient of objective function
        """
        assert exists(path)

        # check parameters
        if 'OPTIMIZE' not in PATH:
            raise ParameterError(PATH, 'OPTIMIZE')

        # check input arguments
        if not exists(path):
            raise Exception()

        self.combine_kernels(path)
        self.process_kernels(path)

        # convert from relative to absolute perturbations
        g = solver.merge(solver.load(
                path +'/'+ 'kernels/sum',
                suffix='_kernel',
                verbose=True))

        g *= solver.merge(solver.load(
                path +'/'+ 'model'))

        # normalize gradient
        if float(PAR.FACTOR) == 1.:
            pass
        elif not PAR.FACTOR:
            pass
        else:
            g *= PAR.FACTOR

        # write gradient
        solver.save(PATH.GRAD +'/'+ 'gradient', 
                    solver.split(g), 
                    suffix='_kernel')

        if PAR.REGULARIZE:
            self.regularize(path)

        if PATH.PRECOND:
            self.precondition(path)

        savenpy(PATH.OPTIMIZE +'/'+ 'g_new',
                solver.merge(solver.load(path +'/'+ 'gradient', suffix='_kernel')))


    def combine_kernels(self, path):
        """ Sums individual source contributions
        """
        system.run('solver', 'combine',
                   hosts='head',
                   path=path +'/'+ 'kernels')


    def process_kernels(self, path):
        """ Performs scaling, smoothing, and preconditioning operations in
            accordance with parameter settings
        """
        # apply clipping
        if PAR.CLIP > 0.:
            system.run('solver', 'clip',
                       hosts='head',
                       path=path + '/' + 'kernels/sum',
                       thresh=PAR.CLIP)

        # apply smoothing
        if PAR.SMOOTH > 0.:
            system.run('solver', 'smooth',
                       hosts='head',
                       path=path + '/' + 'kernels/sum',
                       span=PAR.SMOOTH)


    def precondition(self, path):
        g = solver.merge(solver.load(path +'/'+ 'gradient'), suffix='_kernel')
        g *= solver.merge(solver.load(PATH.PRECOND))

        src = path +'/'+ 'gradient'
        dst = path +'/'+ 'gradient_noprecond'
        unix.mv(src, dst)

        solver.save(path +'/'+ 'gradient',
                    solver.split(g),
                    suffix='_kernel')



    def regularize(self, path):
        raise NotImplementedError
