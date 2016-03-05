
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.array import grid2mesh, mesh2grid, stack
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import
from seisflows.tools.math import nabla


PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class tikhonov0(custom_import('postprocess', 'regularize')):
    """ Adds regularization options to base class

        Available options include 0-, 1-, and 2- order Tikhonov and total
        variation regularization. While the underlying theory is classical,
        application to unstructured numerical grids via the 
        "seisflows.tools.math.nabla" operator is somewhat complicated. 

        So far, can only be used for 2D inversion, because the required spatial
        derivative operator "nabla" is not yet available for 3D grids.
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(tikhonov0, self).check()

        if not PAR.LAMBDA:
            raise ValueError


    def nabla(self, mesh, m, g):
        return m/np.mean(m)

