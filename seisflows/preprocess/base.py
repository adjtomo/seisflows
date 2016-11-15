
import sys
import numpy as np
import obspy

from seisflows.tools import msg, unix
from seisflows.tools.code import exists, Struct
from seisflows.config import   \
    ParameterError

from seisflows.plugins import adjoint, misfit, readers, writers
from seisflows.tools.signal import mute_early, mute_late

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class base(object):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'MISFIT' not in PAR:
            setattr(PAR, 'MISFIT', None)

        if 'BACKPROJECT' not in PAR:
            setattr(PAR, 'BACKPROJECT', None)

        if 'FORMAT' not in PAR:
            raise ParameterError(PAR, 'FORMAT')

        if 'NORMALIZE' not in PAR:
            setattr(PAR, 'NORMALIZE', 'L2')

        if 'MUTE' not in PAR:
            setattr(PAR, 'MUTE', None)

        if 'FILTER' not in PAR:
            setattr(PAR, 'FILTER', None)

        # assertions
        if PAR.FORMAT not in dir(readers):
            print msg.ReaderError
            raise ParameterError()

        if PAR.FORMAT not in dir(writers):
            print msg.WriterError
            raise ParameterError()

        self.check_filter()
        self.check_mute()
        self.check_normalize()


    def setup(self):
        """ Sets up data preprocessing machinery
        """
        # define misfit function and adjoint trace generator
        if PAR.MISFIT:
            self.misfit = getattr(misfit, PAR.MISFIT)
            self.adjoint = getattr(adjoint, PAR.MISFIT)
        elif PAR.BACKPROJECT:
            self.adjoint = getattr(adjoint, PAR.BACKPROJECT)

        # define seismic data reader and writer
        self.reader = getattr(readers, PAR.FORMAT)
        self.writer = getattr(writers, PAR.FORMAT)


    def prepare_eval_grad(self, path='.'):
        """ Prepares solver for gradient evaluation by writing residuals and
          adjoint traces
        """
        solver = sys.modules['seisflows_solver']

        for filename in solver.data_filenames:
            obs = self.reader(path+'/'+'traces/obs', filename)
            syn = self.reader(path+'/'+'traces/syn', filename)

            # process observations
            obs = self.apply_filter(obs)
            obs = self.apply_mute(obs)
            obs = self.apply_normalize(obs)

            # process synthetics
            syn = self.apply_filter(syn)
            syn = self.apply_mute(syn)
            syn = self.apply_normalize(syn)

            if PAR.MISFIT:
                self.write_residuals(path, syn, obs)

            self.write_adjoint_traces(path+'/'+'traces/adj', syn, obs, filename)


    def write_residuals(self, path, syn, dat):
        """ Computes residuals from observations and synthetics
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nn, _ = self.get_network_size(syn)

        filename = path +'/'+ 'residuals'
        if exists(filename):
            rsd = list(np.loadtxt(filename))
        else:
            rsd = []

        for ii in range(nn):
            rsd.append(self.misfit(syn[ii].data, dat[ii].data, nt, dt))

        np.savetxt(filename, rsd)


    def write_adjoint_traces(self, path, syn, dat, channel):
        """ Generates adjoint traces from observed and synthetic traces
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nn, _ = self.get_network_size(syn)

        # NOTE: overwrites syn[:].data
        adj = syn

        for ii in range(nn):
            adj[ii].data = self.adjoint(syn[ii].data, dat[ii].data, nt, dt)

        #self.apply_filter_backwards(adj)

        self.writer(adj, path, channel)


    ### signal processing

    def apply_filter(self, traces):
        if not PAR.FILTER:
            return traces

        elif PAR.FILTER == 'Bandpass':
            traces = _signal.detrend(traces)
            for tr in traces:
                tr.filter('bandpass', freqmin=PAR.FREQMIN, freqmax=PAR.FREQMAX)

                # workaround obspy dtype conversion
                tr.data = tr.data.astype(np.float32)

        elif PAR.FILTER == 'Lowpass':
            traces = _signal.detrend(traces)
            for tr in traces:
                tr.filter('lowpass', freq=PAR.FREQ)

                # workaround obspy dtype conversion
                tr.data = tr.data.astype(np.float32)

        elif PAR.FILTER == 'Highpass':
            traces = _signal.detrend(traces)
            for tr in traces:
                tr.filter('highpass', freq=PAR.FREQ)

                # workaround obspy dtype conversion
                tr.data = tr.data.astype(np.float32)

        else:
            raise ParameterError()

        return traces


    def apply_mute(self, traces):
        if not PAR.MUTE:
            return traces

        elif PAR.MUTE == 'MuteEarlyArrivals':
            # mutes early arrivals
            return mute_early(traces, 
                PAR.MUTESLOPE, # (units: time/distance)
                PAR.MUTECONST, # (units: time)
                self.get_time_scheme(traces),
                self.get_source_coords(traces),
                self.get_receiver_coords(traces))

        elif PAR.MUTE == 'MuteLateArrivals':
            return mute_late(traces,
                PAR.MUTESLOPE, # (units: time/distance)
                PAR.MUTECONST, # (units: time)
                self.get_time_scheme(traces),
                self.get_source_coords(traces),
                self.get_receiver_coords(traces))

        else:
            raise ParameterError()


    def apply_normalize(self, traces):
        if not PAR.NORMALIZE:
            return traces

        elif PAR.NORMALIZE == 'L1':
            # normalize each trace by its L1 norm
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=1)
                if w > 0:
                    tr.data /= w
            return traces

        elif PAR.NORMALIZE == 'L2':
            # normalize each trace by its L2 norm
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=2)
                if w > 0:
                    tr.data /= w
            return traces

        elif PAR.NORMALIZE == 'L2_squared':
            # normalize each trace by its L2 norm squared
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=2)
                if w > 0:
                    tr.data /= w**2.
            return traces

        elif PAR.NORMALIZE == 'L2_summed':
            # normalize all traces by their combined L2 norm
            w = 0.
            for tr in traces:
                w += np.linalg.norm(tr.data, ord=2)**2.
            for tr in traces:
                tr.data /= w**0.5
            return traces

        else:
            raise ParameterError()


    def apply_filter_backwards(self, traces):
        for tr in traces:
            tr.data = np.flip(tr.data)

        traces = self.apply_filter()

        for tr in traces:
            tr.data = np.flip(tr.data)

        return traces



    ### additional parameter checking

    def check_filter(self):
        """ Checks filter settings
        """
        if not PAR.FILTER:
            pass

        elif PAR.FILTER == 'Bandpass':
            if 'FREQMIN' not in PAR: raise ParameterError('FREQMIN')
            if 'FREQMAX' not in PAR: raise ParameterError('FREQMAX')
            assert 0 < PAR.FREQMIN
            assert PAR.FREQMIN < PAR.FREQMIN
            assert PAR.FREQMAX < infinity

        elif PAR.FILTER == 'Lowpass':
            raise NotImplementedError
            if 'FREQ' not in PAR: raise ParameterError('FREQ')
            assert 0 < PAR.FREQ <= infinity

        elif PAR.FILTER == 'Highpass':
            raise NotImplementedError
            if 'FREQ' not in PAR: raise ParameterError('FREQ')
            assert 0 <= PAR.FREQ < infinity

        else:
            raise ParameterError()


    def check_mute(self):
        """ Checks mute settings
        """
        if not PAR.MUTE:
            pass

        elif PAR.MUTE in ['MuteEarlyArrivals', 'MuteLateArrivals']:
            assert 'MUTESLOPE' in PAR
            assert 'MUTECONST' in PAR
            assert PAR.MUTESLOPE >= 0.

        else:
            raise ParameterError()


    def check_normalize(self):
        if not PAR.NORMALIZE:
            pass
        elif PAR.NORMALIZE in ['L1', 'L2', 'L2_squared']:
            pass
        else:
            raise ParameterError()



    ### utility functions

    def get_time_scheme(self, traces):
        # FIXME: extract time scheme from trace headers rather than parameters file
        nt = PAR.NT
        dt = PAR.DT
        t0 = 0.
        return nt, dt, t0


    def get_network_size(self, traces):
        nrec = len(traces)
        nsrc = 1
        return nrec, nsrc


    def get_receiver_coords(self, traces):
        if PAR.FORMAT in ['SU', 'su']:
            rx = []
            ry = []
            rz = []
            for trace in traces:
                rx += [trace.stats.su.trace_header.group_coordinate_x]
                ry += [trace.stats.su.trace_header.group_coordinate_y]
                rz += [0.]
            return rx, ry, rz

        else:
             raise NotImplementedError


    def get_source_coords(self, traces):
        if PAR.FORMAT in ['SU', 'su']:
            sx = []
            sy = []
            sz = []
            for trace in traces:
                sx += [trace.stats.su.trace_header.source_coordinate_x]
                sy += [trace.stats.su.trace_header.source_coordinate_y]
                sz += [0.]
            return sx, sy, sz

        else:
             raise NotImplementedError


