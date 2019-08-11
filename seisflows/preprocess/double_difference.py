#
# This is Seisflows
#
# See LICENCE file
#
###############################################################################

# Import system modules
import sys
from os.path import exists

# Import Numpy and Obspy
import numpy as np
from obspy.core import Stream, Trace

# Local imports
from seisflows.plugins import adjoint, misfit
from seisflows.tools import unix
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']


class double_difference(custom_import('preprocess', 'base')):
    """ Double-difference data processing class

      Adds double-difference data misfit functions to base class
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(double_difference, self).check()

        if not hasattr(PAR, 'DISTMAX'):
            setattr(PAR, 'DISTMAX', float("inf"))

        if not hasattr(PAR, 'UNITS'):
            setattr(PAR, 'UNITS', 'lonlat')

        if not hasattr(PATH, 'WEIGHTS'):
            setattr(PATH, 'WEIGHTS', None)

        if PATH.WEIGHTS:
            assert exists(PATH.WEIGHTS)

        assert PAR.MISFIT in [
            'Traveltime',
            'TraveltimeInexact']

    def write_residuals(self, path, syn, dat):
        """ Computes residuals from observations and synthetics
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nr, _ = self.get_network_size(syn)
        rx, ry, rz = self.get_receiver_coords(syn)

        dist = np.zeros((nr, nr))
        count = np.zeros(nr)
        delta_syn = np.zeros((nr, nr))
        delta_obs = np.zeros((nr, nr))

        # calculate distances between stations
        for i in range(nr):
            for j in range(i):
                dist[i, j] = self.distance(rx[i], rx[j], ry[i], ry[j])

        # calculate traveltime lags between stations pairs
        for i in range(nr):
            for j in range(i):
                if dist[i, j] > PAR.DISTMAX:
                    continue
                delta_syn[i, j] = self.misfit(syn[i].data, syn[j].data, nt, dt)
                delta_obs[i, j] = self.misfit(dat[i].data, dat[j].data, nt, dt)
                delta_syn[j, i] = -delta_syn[i, j]
                delta_obs[j, i] = -delta_obs[i, j]
                count[i] += 1

        np.savetxt(path + '/' + 'dist_ij', dist)
        np.savetxt(path + '/' + 'count', count)
        np.savetxt(path + '/' + 'delta_syn_ij', delta_syn)
        np.savetxt(path + '/' + 'delta_obs_ij', delta_obs)
        np.savetxt(path + '/' + 'rsd_ij', delta_syn-delta_obs)

        # to get residuals, sum over all station pairs
        rsd = abs(delta_syn-delta_obs).sum(axis=0)

        # apply optional weights
        if PATH.WEIGHTS:
            rsd *= self.load_weights()

        # write residuals to text file
        filename = path + '/' + 'residuals'
        if exists(filename):
            rsdlist = list(np.loadtxt(filename))
        else:
            rsdlist = []
        rsdlist += [rsd]
        np.savetxt(filename, rsdlist)

    def sum_residuals(self):
        """ Sums squares of residuals
        """
        total_misfit = 0.
        for path in paths:
            total_misfit += np.sum(np.loadtxt(path)**2.)
        return total_misfit

    def write_adjoint_traces(self, path, syn, dat, channel):
        """ Computes adjoint traces from observed and synthetic traces
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nr, _ = self.get_network_size(syn)

        Del = np.loadtxt(path + '/' + '../../delta_syn_ij')
        rsd = np.loadtxt(path + '/' + '../../rsd_ij')

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

                adj[i].data += rsd[i, j] * \
                               self.adjoint_dd(si, sj, +Del[i, j], nt, dt)
                adj[j].data -= rsd[i, j] * \
                               self.adjoint_dd(sj, si, -Del[i, j], nt, dt)

        # optional weighting
        adj = self.apply_weights(adj)

        # write adjoint traces
        self.writer(adj, path, channel)

    def adjoint_dd(self, si, sj, t0, nt, dt):
        """ Returns contribution to adjoint source from a single double
            difference measurement
        """
        vi = np.zeros(nt)
        vj = np.zeros(nt)

        vi[1:-1] = (si[2:] - si[0:-2])/(2.*dt)
        vj[1:-1] = (sj[2:] - sj[0:-2])/(2.*dt)

        vjo = self.shift(vj, -t0/dt)

        w = sum(vi*vjo*dt)
        w = max(vjo)
        if w:
            vjo /= w

        return vjo

    def apply_weights(self, traces):
        if not PATH.WEIGHTS:
            return traces

        else:
            w = self.load_weights()
            for i, trace in enumerate(traces):
                trace.data *= w[i]
            return traces

    def load_weights(self):
            return np.loadtxt(PATH.WEIGHTS)[:, -1]

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

    def distance(self, x1, y1, x2, y2):
        if PAR.UNITS in ['lonlat']:
            dlat = np.radians(y2-y1)
            dlon = np.radians(x2-x1)
            a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(y1)) \
                * np.cos(np.radians(y2)) * np.sin(dlon/2) * np.sin(dlon/2)
            D = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
            D *= 180/np.pi
            return D

        else:
            return ((x1-x2)**2 + (y1-y2)**2)**0.5
