
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


class regularize(loadclass('postprocess', 'base')):
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
        super(regularize, self).check()

        if 'RADIUS' not in PAR:
            raise ParameterError(PAR, 'RADIUS')

        if 'LAMBDA' not in PAR:
            setattr(PAR, 'LAMBDA', 0.)


    def setup(self):
        """ Performs any required initialization or setup tasks
        """
        pass


    def process_kernels(self, path):
        """ Masks source and receiver artifacts
        """
        fullpath = path +'/'+ 'kernels'
        assert exists(path)

        if exists(fullpath +'/'+ 'sum'):
            unix.mv(fullpath +'/'+ 'sum', fullpath +'/'+ 'sum_nomask')

        # mask sources and receivers
        system.run('postprocess', 'mask', 
                   hosts='all', 
                   path=fullpath)

        system.run('solver', 'combine',
                   hosts='head',
                   path=fullpath)


    def mask(self, path=''):
        import preprocess
        preprocess.setup()

        fullpath = path +'/'+  solver.getname
        g = solver.load(fullpath, suffix='_kernel')
        if not PAR.RADIUS:
            return

        try:
            x = g['x'][0]
            z = g['z'][0]
            mesh = stack(x, z)
        except:
            from seisflows.seistools.io import loadbin
            model_path = PATH.OUTPUT +'/'+ 'model_true'
            x = loadbin(model_path, 0, 'x')
            z = loadbin(model_path, 0, 'z')
            mesh = stack(x, z)

        lx = x.max() - x.min()
        lz = z.max() - z.min()
        nn = x.size
        nx = np.around(np.sqrt(nn*lx/lz))
        nz = np.around(np.sqrt(nn*lz/lx))
        dx = lx/nx
        dz = lz/nz

        sigma = 0.5*PAR.RADIUS*(dx+dz)
        _, h = preprocess.load(solver.getpath +'/'+ 'traces/obs')

        # mask sources
        mask = np.exp(-0.5*((x-h.sx[0])**2.+(z-h.sy[0])**2.)/sigma**2.)
        for key in solver.parameters:
            weight = np.sum(mask*g[key][0])/np.sum(mask)
            g[key][0] *= 1.-mask
            g[key][0] += mask*weight

        # mask receivers
        for ir in range(h.nr):
            mask = np.exp(-0.5*((x-h.rx[ir])**2.+(z-h.ry[ir])**2.)/sigma**2.)
            for key in solver.parameters:
                weight = np.sum(mask*g[key][0])/np.sum(mask)
                g[key][0] *= 1.-mask
                g[key][0] += mask*weight

        solver.save(fullpath, g, suffix='_kernel')


    def regularize(self, path):
        assert (exists(path))

        g = solver.load(path +'/'+ 'gradient', suffix='_kernel')
        if not PAR.LAMBDA:
            return solver.merge(g)

        m = solver.load(path +'/'+ 'model')
        try:
            x = m['x'][0]
            z = m['z'][0]
            mesh = stack(x, z)
        except:
            from seisflows.seistools.io import loadbin
            model_path = PATH.OUTPUT +'/'+ 'model_true'
            x = loadbin(model_path, 0, 'x')
            z = loadbin(model_path, 0, 'z')
            mesh = stack(x, z)

        for key in solver.parameters:            
            for iproc in range(PAR.NPROC):
                g[key][iproc] += PAR.LAMBDA *\
                    self.nabla(mesh, m[key][iproc], g[key][iproc])

        return solver.merge(g)


    def nabla(self, mesh, m, g):
        raise NotImplementedError("Must be implemented by subclass.")

