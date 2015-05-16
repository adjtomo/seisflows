
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import postprocess


class test_postprocess(object):
    """ Postprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'KERNELS' not in PATH:
            raise ParameterError(PATH, 'KERNELS')

        if not exists(PATH.KERNELS +'/'+ 'kernels'):
            raise Exception('Bad path: ' % PATH.KERNELS)


        # check postprocessing settings
        if 'SCALE' not in PAR:
            setattr(PAR, 'SCALE', False)

        if 'CLIP' not in PAR:
            setattr(PAR, 'CLIP', 0.)

        if 'SMOOTH' not in PAR:
            setattr(PAR, 'SMOOTH', 0.)

        if 'PRECOND' not in PATH:
            setattr(PATH, 'PRECOND', None)


    def main(self):
        """ Writes gradient of objective function
        """
        #print 'Combining kernels...'
        #postprocess.combine_kernels(PATH.KERNELS)

        print 'Processing kernels...'
        postprocess.process_kernels(PATH.KERNELS)

