
import sys
import numpy as np

from os.path import exists
from obspy.core import Stream, Trace

from seisflows.plugins import adjoint, misfit
from seisflows.tools import unix
from seisflows.tools.tools import Struct
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']


class double_difference(custom_import('preprocess', 'base')):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(double_difference, self).check()

        if PAR.MISFIT not in ['Traveltime']:
            raise Exception

        if not hasattr(PAR, 'DEBUG_INEXACT_CC'):
            setattr(PAR, 'DEBUG_INEXACT_CC', False)


    def write_residuals(self, path, syn, dat):
        """ Computes residuals from observations and synthetics
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nr, _ = self.get_network_size(syn)

        # compute traveltime differences between stations
        d = np.zeros((nr, nr))
        dd = np.zeros((nr, nr))

        for i in range(nr):
            for j in range(i):

                if PAR.VERBOSE >= 2:
                    if system.getnode() == 0:
                        print i,j

                d[i,j] = self.misfit_dd(syn[i].data, syn[j].data, nt, dt)
                dd[i,j] = d[i,j] - self.misfit_dd(dat[i].data, dat[j].data, nt, dt)

        np.savetxt(path +'/'+ 'dij', d)
        np.savetxt(path +'/'+ 'ddij', dd)

        filename = path +'/'+ 'residuals'
        if exists(filename):
            rsd = list(np.loadtxt(filename))
        else:
            rsd = []

        # to get residual, sum over all station pairs
        for i in range(nr):
            rsd += [np.sum(abs(dd), 0)]
        np.savetxt(filename, rsd)


    def write_adjoint_traces(self, path, syn, dat, channel):
        """ Computes adjoint traces from observed and synthetic traces
        """

        nt, dt, _ = self.get_time_scheme(syn)
        nr, _ = self.get_network_size(syn)

        d = np.loadtxt(path +'/'+ '../../dij')
        dd = np.loadtxt(path +'/'+ '../../ddij')

        # initialize adjoint traces
        adj = Stream()
        for i in range(nr):
            adj.append(Trace(
                data=np.zeros(nt, dtype='float32'),
                header=syn[i].stats))

        # generate adjoint traces
        for i in range(nr):
            for j in range(i):

                if PAR.VERBOSE > 2:
                    print i,j

                ai = adj[i].data
                aj = adj[j].data
                si = syn[i].data
                sj = syn[j].data

                ai_tmp = self.adjoint_dd(si, sj, +d[i,j], nt, dt)
                aj_tmp = self.adjoint_dd(sj, si, -d[i,j], nt, dt)

                ai -= dd[i,j] * ai_tmp
                aj += dd[i,j] * aj_tmp

        # write adjoint traces
        self.writer(adj, path, channel)


    def adjoint_dd(self, si, sj, t0, nt, dt):
        vi = np.zeros(nt)
        vj = np.zeros(nt)

        vi[1:-1] = (si[2:] - si[0:-2])/(2.*dt)
        vj[1:-1] = (sj[2:] - sj[0:-2])/(2.*dt)

        vjo = self.shift(vj, t0/dt)

        w  = sum(vi*vjo*dt)
        w = max(vjo)
        if w:
            vjo /= w

        return vjo


    def misfit_dd(self, si, sj, nt, dt):
        if PAR.DEBUG_INEXACT_CC:
            # much faster but possibly inaccurate
            itmax = np.argmax(si)
            jtmax = np.argmax(sj)
            return (itmax-jtmax)*dt
        else:
            cc = np.convolve(si, np.flipud(sj))
            it = np.argmax(cc)
            return (it-nt+1)*dt


    def shift(self, v, it):
        if it == 0:
            return v

        nt = len(v)
        vo = np.zeros(nt)
        if it > 0:
            # shift right
            vo[it:] = v[:-it]
        else:
            # shift left
            vo[:it] = v[-it:]
        return vo

            
