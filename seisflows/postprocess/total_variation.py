
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.array import grid2mesh, mesh2grid, stack
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import
from seisflows.tools.math import nabla, tv


PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class total_variation(custom_import('postprocess', 'regularize')):
    """ Adds regularization options to base class

        So far, can only be used for 2D inversion, because the required spatial
        derivative operator "nabla" is not yet available for 3D grids.
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(total_variation, self).check()

        if not PAR.LAMBDA:
            raise ValueError

        if not hasattr(PAR, 'EPSILON'):
            setattr(PAR, 'EPSILON', 0.)


    def nabla(self, mesh, m, g):
        M, grid = mesh2grid(g, mesh)
        DM = tv(M, epsilon=PAR.EPSILON)
        dm = grid2mesh(DM, grid, mesh)
        return dm/np.mean(m)

