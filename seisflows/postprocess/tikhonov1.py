
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.array import grid2mesh, mesh2grid, stack
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass
from seisflows.tools.math import nabla


PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class tikhonov1(loadclass('postprocess', 'regularize')):
    """ Adds regularization options to base class

        Available options include 0-, 1-, and 2- order Tikhonov and total
        variation regularization. While the underlying theory is classical,
        these options are experimental in the sense that their application to
        unstructured numerical grids is quite new.

        SO FAR, CAN ONLY BE USED FOR 2D WAVEFORM INVERSION.
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(tikhonov1, self).check()

        if 'UPDATETYPE' not in PAR:
            setattr(PAR, 'UPDATETYPE', 'Creeping')

        #if not PAR.LAMBDA:
        #    raise ValueError


    def nabla(self, mesh, m, g):
        if PAR.UPDATETYPE in ['Creeping']:
            G, grid = mesh2grid(g, mesh)
            DG = nabla(G, order=1)
            dg = grid2mesh(DG, grid, mesh)
            return -dg/np.mean(m)

        elif PAR.UPDATETYPE in ['Jumping']:
            M, grid = mesh2grid(m, mesh)
            DM = nabla(M, order=1)
            dm = grid2mesh(DM, grid, mesh)
            return dm/np.mean(m)

