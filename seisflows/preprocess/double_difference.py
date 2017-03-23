
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

        # specify defaults 
        if not hasattr(PAR, 'DISTMAX'):
            setattr(PAR, 'DISTMAX', np.inf)

        if not hasattr(PAR, 'WEIGHTS'):
            setattr(PAR, 'WEIGHTS', None)

        # error checking
        if PAR.MISFIT not in ['Traveltime', 'TraveltimeInexact']:
            raise Exception


    def write_residuals(self, path, syn, dat):
        """ Computes residuals from observations and synthetics
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nr, _ = self.get_network_size(syn)
        rx, ry, rz = self.get_receiver_coords(syn)

        # calculate distances between stations
        dist = np.zeros((nr,nr))
        count = np.zers(nr)
        for i in range(nr):
            for j in range(i):
                dist[i,j] = ((rx[i]-rx[j])**2 +
                             (ry[i]-ry[j])**2 +
                             (rz[i]-rz[j])**2)**0.5

        # calculate traveltime differences between stations
        delta_syn = np.zeros((nr,nr))
        delta_obs = np.zeros((nr,nr))

        for i in range(nr):
            for j in range(i):
                if dist[i,j] > PAR.DISTMAX: 
                    continue

                delta_syn[i,j] = self.get_misfit(syn[i].data,syn[j].data,nt,dt)
                delta_obs[i,j] = self.get_misfit(dat[i].data,dat[j].data,nt,dt)
                count[i] += 1

        np.savetxt(path +'/'+ 'dist_ij', dist)
        np.savetxt(path +'/'+ 'count', count)
        np.savetxt(path +'/'+ 'delta_syn_ij', delta_syn)
        np.savetxt(path +'/'+ 'delta_obs_ij', delta_obs)
        np.savetxt(path +'/'+ 'rsd_ij', delta_syn-delta_obs)

        # to get residual, sum over all station pairs
        filename = path +'/'+ 'residuals'
        if exists(filename):
            rsd = list(np.loadtxt(filename))
        else:
            rsd = []
        rsd += [np.sum(abs(delta_syn-delta_obs), 0)]
        np.savetxt(filename, rsd)


    def write_adjoint_traces(self, path, syn, dat, channel):
        """ Computes adjoint traces from observed and synthetic traces
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nr, _ = self.get_network_size(syn)

        Del = np.loadtxt(path +'/'+ '../../delta_syn_ij')
        rsd = np.loadtxt(path +'/'+ '../../rsd_ij')

        # initialize trace arrays
        adj = Stream()
        for i in range(nr):
            adj.append(Trace(
                data=np.zeros(nt, dtype='float32'),
                header=syn[i].stats))

        # generate adjoint traces
        for i in range(nr):
            for j in range(i):
                si = syn[i].data
                sj = syn[j].data

                adj[i].data -= rsd[i,j] * \
                               self.get_adjoint(si, sj, +Del[i,j], nt, dt)
                adj[j].data += rsd[i,j] * \
                               self.get_adjoint(sj, si, -Del[i,j], nt, dt)

        # apply weights
        adj = self.apply_weights(adj)

        # write adjoint traces
        self.writer(adj, path, channel)


    def get_adjoint(self, si, sj, t0, nt, dt):
        """ Returns contribution to adjoint source from a single double 
         difference measurement
        """
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


    def get_misfit(self, si, sj, nt, dt):
        """ Calculates misfit for a single double difference measurement
        """
        if PAR.MISFIT == 'TraveltimeInexact':
            # much faster but possibly inaccurate
            itmax = np.argmax(si)
            jtmax = np.argmax(sj)
            return (itmax-jtmax)*dt

       elif PAR.MISFIT == 'Traveltime'
            cc = np.convolve(si, np.flipud(sj))
            it = np.argmax(cc)
            return (it-nt+1)*dt

       else:
           raise ParameterError


    def apply_weights(self, traces):
        if not PAR.WEIGHTS:
            return traces

        elif PAR.WEIGHTS in ['Simple']:
            w = np.loadtxt('count')
            for i,trace in enumerate(traces):
                trace.data /= w[i]

        elif PAR.WEIGHTS in ['UserSupplied']:
            w = np.loadtxt(PATH.WEIGHTS)
            for i,trace in enumerate(traces):
                trace.data /= w[i]

        elif PAR.WEIGHTS in ['Geographic']:
            raise NotImplementedError


    def shift(self, v, it):
        """ Shifts time series a given number of steps
        """
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


