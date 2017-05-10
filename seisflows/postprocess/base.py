
import sys
import numpy as np

from os.path import join
from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import ParameterError

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']


class base(object):
    """ Gradient postprocessing class
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

        if 'PRECOND' not in PAR:
            setattr(PAR, 'PRECOND', False)

        # check paths
        if 'MASK' not in PATH:
            setattr(PATH, 'MASK', None)

        if 'PRECOND' not in PATH:
            setattr(PATH, 'PRECOND', None)

        if PATH.MASK:
            assert exists(PATH.MASK)


    def setup(self):
        """ Performs any required initialization or setup tasks
        """
        pass


    def write_gradient(self, path):
        """ Reads kernels and writes gradient of objective function
        """
        if not exists(path):
            raise Exception()

        system.run('postprocess', 'process_kernels',
                 hosts='head',
                 path=path,
                 parameters=solver.parameters)

        g = solver.merge(solver.load(
                 path +'/'+ 'kernels/sum',
                 suffix='_kernel'))

        if PAR.LOGARITHMIC:
            # convert from logarithmic to absolute perturbations
            g *= solver.merge(solver.load(path +'/'+ 'model'))
        self.save(path, g)

        if PATH.MASK:
            # apply mask
            g *= solver.merge(solver.load(PATH.MASK))
            self.save(path, g, backup='nomask')


    def process_kernels(self, path='', parameters=[]):
        solver.combine(
               path +'/'+ 'kernels',
               parameters)

        if PAR.SMOOTH > 0.:
            solver.smooth(
                   path=path +'/'+ 'kernels/sum',
                   parameters=parameters,
                   span=PAR.SMOOTH)


    def save(self, path, g, backup=None):
        if backup:
            src = path +'/'+ 'gradient'
            dst = path +'/'+ 'gradient_'+backup
            unix.mv(src, dst)

        solver.save(solver.split(g),
                    path +'/'+ 'gradient',
                    parameters=solver.parameters,
                    suffix='_kernel')


