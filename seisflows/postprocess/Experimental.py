
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.array import grid2mesh, mesh2grid, nabla, stack
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ParameterObj, ParameterError

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import solver


class Experimental(loadclass('postprocess', 'base')):
    """ Adds alternate regularization options to base class

        Regularization strategies are experimental in the sense that 
        we are trying to make them perform better when the model is 
        expressed on an unstructured numerical grid.
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(Experimental, self).check()

        if 'REGULARIZE' not in PAR:
            setattr(PAR, 'REGULARIZE', None)

        if 'LAMBDA' not in PAR:
            setattr(PAR, 'LAMBDA', 0.)


    def process_kernels(self, path):
        """ Performs scaling, smoothing, and preconditioning operations in
            accordance with parameter settings
        """
        assert (exists(path))

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

        # regularize
        if PAR.REGULARIZE:
            parts = solver.load(path +'/'+ 'kernels/sum')
            assert 'x' in parts
            assert 'z' in parts
            x = parts['x'][0]
            z = parts['z'][0]
            mesh = stack(x, z)

            for key in solver.parameters:            
                for iproc in range(PAR.NPROC):
                    parts[key][iproc] += PAR.LAMBDA*self.regularize(parts[key][0], mesh)

            src = path +'/'+ 'kernels/sum'
            dst = path +'/'+ 'kernels/sum_noregularize'
            unix.mv(src, dst)

            solver.save(path +'/'+ 'kernels/sum',
                        parts,
                        suffix='_kernel')
 
        # precondition
        if PATH.PRECOND:
            g = solver.merge(solver.load(path +'/'+ 'kernels/sum'))
            g *= solver.merge(solver.load(PATH.PRECOND))

            src = path +'/'+ 'kernels/sum'
            dst = path +'/'+ 'kernels/sum_noprecond'
            unix.mv(src, dst)

            solver.save(path +'/'+ 'kernels/sum',
                        solver.split(g),
                        suffix='_kernel')


    def regularize(self, v, mesh):
        if PAR.REGULARIZE in ['TotalVariation']:
            V, grid = mesh2grid(v, mesh)
            print V.min(), V.max()
            DV = nabla(V, order=1)
            print DV.min(), DV.max()
            dv = grid2mesh(DV, grid, mesh)
            print dv.min(), dv.max()
            return dv

        elif PAR.REGULARIZE in ['Tikhonov0', 'Damping']:
            return v

        elif PAR.REGULARIZE in ['Tikhonov1', 'Smoothing']:
            V, grid = mesh2grid(v, mesh)
            DV = nabla(V, order=1)
            dv = grid2mesh(DV, grid, mesh)
            return dv

        elif PAR.REGULARIZE in ['Tikhonov2']:
            V, grid = mesh2grid(v, mesh)
            DV = nabla(V, order=2)
            dv = grid2mesh(DV, grid, mesh)
            return dv



    def test_regularize(self, model):
        # 'Total Variation'
            V, grid = mesh2grid(v, mesh)
            DV = nabla(V, order=1)
            dv = grid2mesh(DV, grid)

        # 'Tikhonov1', 'Smoothing'
            V, grid = mesh2grid(v, mesh)
            DV = nabla(V, order=1)
            dv = grid2mesh(DV, grid)
            return dv

