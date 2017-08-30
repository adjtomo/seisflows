
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

        if 'KERNELTYPE' not in PAR:
            setattr(PAR, 'KERNELTYPE', 'Relative')

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
        """ Can be used to perform any required initialization or setup tasks
        """
        pass


    def write_gradient(self, path):
        """ Processes and combines contributions to the gradient from
          individual sources
        """
        if not exists(path):
            raise Exception

        system.run('postprocess', 'process_kernels',
                 hosts='head',
                 path=path+'/kernels',
                 parameters=solver.parameters)

        g = solver.merge(solver.load(
                 path +'/'+ 'kernels/sum',
                 suffix='_kernel'))

        self.save(g, path)

        if PAR.KERNELTYPE=='Relative':
            # convert from relative to absolute perturbations
            g *= solver.merge(solver.load(path +'/'+ 'model'))
            self.save(g, path, backup='relative')

        if PATH.MASK:
            # apply mask
            g *= solver.merge(solver.load(PATH.MASK))
            self.save(g, path, backup='nomask')


    def process_kernels(self, path='', parameters=[]):
        """ Combines contributions from individual sources and performs any 
         required processing steps

          INPUT
              PATH - directory containing sensitivity kernels
              PARAMETERS - list of material parameters to be operated on
        """
        if not exists(path):
            raise Exception

        if not parameters:
            parameters = solver.parameters

        solver.combine(
               input_path=path,
               output_path=path+'/'+'sum',
               parameters=parameters)

        if PAR.SMOOTH > 0.:
            src = path+'/'+'sum'
            dst = path+'/'+'sum_nosmooth' 
            unix.mv(src, dst)

            solver.smooth(
                   input_path=path+'/'+'sum_nosmooth',
                   output_path=path+'/'+'sum',
                   parameters=parameters,
                   span=PAR.SMOOTH)


    def save(self, g, path='', parameters=[], backup=None):
        """ Convience function for saving dictionary representation of 
          gradient
        """
        if not exists(path):
            raise Exception

        if not parameters:
            parameters = solver.parameters

        if backup:
            src = path +'/'+ 'gradient'
            dst = path +'/'+ 'gradient_'+backup
            unix.mv(src, dst)

        solver.save(solver.split(g),
                    path +'/'+ 'gradient',
                    parameters=parameters,
                    suffix='_kernel')


