
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
        # check parameters
        if 'CLIP' not in PAR:
            setattr(PAR, 'CLIP', 0.)

        if 'SMOOTH' not in PAR:
            setattr(PAR, 'SMOOTH', 0.)

        if 'LOGARITHMIC' not in PAR:
            setattr(PAR, 'LOGARITHMIC', True)

        if 'REGULARIZE' not in PAR:
            setattr(PAR, 'REGULARIZE', False)

        if 'PRECOND' not in PAR:
            setattr(PAR, 'PRECOND', False)

        # check paths
        if 'MASK' not in PATH:
            setattr(PATH, 'MASK', None)

        if 'PRECOND' not in PATH:
            setattr(PATH, 'PRECOND', None)

        # assertions
        if PAR.PRECOND:
            assert exists(PATH.PRECOND)

        if PATH.PRECOND:
            assert PAR.PRECOND
            assert exists(PATH.PRECOND)

        if PATH.MASK:
            assert exists(PATH.MASK)


    def setup(self):
        """ Performs any required initialization or setup tasks
        """
        pass


    def write_gradient(self, path):
        """ Reads kernels and writes gradient of objective function
        """
        if 'OPTIMIZE' not in PATH:
            raise ParameterError(PATH, 'OPTIMIZE')

        if not exists(path):
            raise Exception()

        self.combine_kernels(path)
        self.process_kernels(path)

        g = solver.merge(solver.load(
                         path +'/'+ 'kernels/sum',
                         suffix='_kernel',
                         verbose=True))

        if PAR.LOGARITHMIC:
            # convert from logarithmic to absolute perturbations
            g *= solver.merge(solver.load(path +'/'+ 'model'))
            self.save(path, g)

        if PATH.MASK:
            # apply mask
            g *= solver.merge(solver.load(PATH.MASK))
            self.save(path, g, backup='nomask')

        if PAR.REGULARIZE:
            # apply regularization
            g = self.apply_regularization(path)
            self.save(path, g, backup='noregularize')

        savenpy(PATH.OPTIMIZE +'/'+ 'g_new', g)


    def write_preconditioner(self):
        """ Reads and writes user-supplied diagonal preconditioner
        """
        if 'OPTIMIZE' not in PATH:
            raise ParameterError(PATH, 'OPTIMIZE')

        savenpy(PATH.OPTIMIZE +'/'+ 'precond', 
                solver.merge(solver.load(PATH.PRECOND)))


    def combine_kernels(self, path):
        system.run('solver', 'combine',
                   hosts='head',
                   path=path +'/'+ 'kernels')


    def process_kernels(self, path):
        if PAR.SMOOTH > 0.:
            system.run('solver', 'smooth',
                       hosts='head',
                       path=path + '/' + 'kernels/sum',
                       span=PAR.SMOOTH)


    def apply_regularization(self, *args, **kwargs):
        raise NotImplementedError("Must be implemented by subclass")


    def save(self, path, v, backup=None):
        if backup:
            src = path +'/'+ 'gradient'
            dst = path +'/'+ 'gradient_'+backup
            unix.mv(src, dst)

        solver.save(path +'/'+ 'gradient',
                    solver.split(v),
                    suffix='_kernel')

