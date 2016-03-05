
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


class regularize(custom_import('postprocess', 'base')):
    """ Adds regularization options to base class

        This parent class is only an abstract base class; see child classes
        TIKHONOV1, TIKHONOV1, and TOTAL_VARIATION for usable regularization.

        Prior to regularizing gradient, near field artifacts must be corrected.
        The "FIXRADIUS" parameter specifies the radius, in number of GLL points,
        within which the correction is applied.
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(regularize, self).check()

        if 'FIXRADIUS' not in PAR:
            setattr(PAR, 'FIXRADIUS', 7.5)

        if 'LAMBDA' not in PAR:
            setattr(PAR, 'LAMBDA', 0.)


    def write_gradient(self, path):
        super(regularize, self).write_gradient(path)

        g = self.regularize(path)
        self.save(path, g, backup='noregularize')


    def process_kernels(self, path, parameters):
        """ Processes kernels in accordance with parameter settings
        """
        fullpath = path +'/'+ 'kernels'
        assert exists(path)

        if exists(fullpath +'/'+ 'sum'):
            unix.mv(fullpath +'/'+ 'sum', fullpath +'/'+ 'sum_nofix')

        # mask sources and receivers
        system.run('postprocess', 'fix_near_field', 
                   hosts='all', 
                   path=fullpath)

        system.run('solver', 'combine',
                   hosts='head',
                   path=fullpath,
                   parameters=parameters)


    def fix_near_field(self, path=''):
        """
        """
        import preprocess
        preprocess.setup()

        name = solver.check_source_names()[solver.getnode]
        fullpath = path +'/'+ name
        g = solver.load(fullpath, suffix='_kernel')
        if not PAR.FIXRADIUS:
            return

        mesh = self.getmesh()
        x,z = self.getxz()

        lx = x.max() - x.min()
        lz = z.max() - z.min()
        nn = x.size
        nx = np.around(np.sqrt(nn*lx/lz))
        nz = np.around(np.sqrt(nn*lz/lx))
        dx = lx/nx
        dz = lz/nz

        sigma = 0.5*PAR.FIXRADIUS*(dx+dz)
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
        mesh = self.getmesh()

        for key in solver.parameters:            
            for iproc in range(PAR.NPROC):
                g[key][iproc] += PAR.LAMBDA *\
                    self.nabla(mesh, m[key][iproc], g[key][iproc])
                    #self.nabla(m[key][iproc], g[key][iproc] , mesh, h)

        return solver.merge(g)


    def nabla(self, mesh, m, g):
        raise NotImplementedError("Must be implemented by subclass.")


    def getmesh(self):
        model_path = PATH.OUTPUT +'/'+ 'model_init'
        try:
            m = solver.load(model_path)
            x = m['x'][0]
            z = m['z'][0]
            mesh = stack(x, z)
        except:
            from seisflows.seistools.io import loadbin
            x = loadbin(model_path, 0, 'x')
            z = loadbin(model_path, 0, 'z')
            mesh = stack(x, z)
        return mesh


    def getxz(self):
        model_path = PATH.OUTPUT +'/'+ 'model_init'
        try:
            m = solver.load(model_path)
            x = m['x'][0]
            z = m['z'][0]
        except:
            from seisflows.seistools.io import loadbin
            x = loadbin(model_path, 0, 'x')
            z = loadbin(model_path, 0, 'z')
        return x,z


      #if PAR.SOLVER in ['specfem2d', 'elastic2d']:
      #    from seisflows.seistools.io import loadbin
      #    x = loadbin(model_path, 0, 'x')
      #    z = loadbin(model_path, 0, 'z')
      #    mesh = stack(x, z)
      #elif PAR.SOLVER in ['specfem2d_legacy']:
      #    m = solver.load(model_path)
      #    x = g['x'][0]
      #    z = g['z'][0]
      #    mesh = stack(x, z)
      #else:
      #    msg = "postprocess.regularize can only be used with solver.specfem2d"
      #    raise Exception(msg)
