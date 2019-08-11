#
# This is Seisflows
#
# See LICENCE file
#
###############################################################################

# Import system modules
import sys

# Import Numpy and Obspy
import numpy as np
import obspy

# Local imports
from seisflows.tools import msg, unix
from seisflows.tools.tools import exists, getset
from seisflows.config import ParameterError
from seisflows.plugins import adjoint, misfit, readers, writers
from seisflows.tools import signal

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class base(object):
    """ Data preprocessing class

      Provides data processing functions for seismic traces, with options for
      data misfit, filtering, normalization and muting
    """

    def check(self):
        """ Checks parameters and paths
        """
        # used for inversion
        if 'MISFIT' not in PAR:
            setattr(PAR, 'MISFIT', None)

        # used for migration
        if 'BACKPROJECT' not in PAR:
            setattr(PAR, 'BACKPROJECT', None)

        # data file format
        if 'FORMAT' not in PAR:
            raise ParameterError(PAR, 'FORMAT')

        # data normalization option
        if 'NORMALIZE' not in PAR:
            setattr(PAR, 'NORMALIZE', None)

        # data muting option
        if 'MUTE' not in PAR:
            setattr(PAR, 'MUTE', None)

        # data filtering option
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
        """
         Prepares solver for gradient evaluation by writing residuals and
         adjoint traces

         :input path: directory containing observed and synthetic seismic data
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

            self.write_adjoint_traces(path+'/'+'traces/adj', syn, obs,
                                      filename)

    def write_residuals(self, path, syn, obs):
        """
        Computes residuals

        :input path: location "adjoint traces" will be written
        :input syn: obspy Stream object containing synthetic data
        :input obs: obspy Stream object containing observed data
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nn, _ = self.get_network_size(syn)

        residuals = []
        for ii in range(nn):
            residuals.append(self.misfit(syn[ii].data, obs[ii].data, nt, dt))

        filename = path+'/'+'residuals'
        if exists(filename):
            residuals.extend(list(np.loadtxt(filename)))

        np.savetxt(filename, residuals)

    def sum_residuals(self, files):
        """
        Sums squares of residuals

        :input files: list of single-column text files containing residuals
        :output total_misfit: sum of squares of residuals
        """
        total_misfit = 0.
        for filename in files:
            total_misfit += np.sum(np.loadtxt(filename)**2.)
        return total_misfit

    def write_adjoint_traces(self, path, syn, obs, channel):
        """
        Writes "adjoint traces" required for gradient computation

        :input path: location "adjoint traces" will be written
        :input syn: obspy Stream object containing synthetic data
        :input obs: obspy Stream object containing observed data
        :input channel: channel or component code used by writer
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nn, _ = self.get_network_size(syn)

        adj = syn
        for ii in range(nn):
            adj[ii].data = self.adjoint(syn[ii].data, obs[ii].data, nt, dt)

        self.writer(adj, path, channel)

    # Signal processing
    def apply_filter(self, traces):
        if not PAR.FILTER:
            return traces

        elif PAR.FILTER == 'Bandpass':
            for tr in traces:
                tr.detrend('demean')
                tr.detrend('linear')
                tr.taper(0.05, type='hann')
                tr.filter('bandpass',
                          zerophase=True,
                          freqmin=PAR.FREQMIN,
                          freqmax=PAR.FREQMAX)

        elif PAR.FILTER == 'Lowpass':
            for tr in traces:
                tr.detrend('demean')
                tr.detrend('linear')
                tr.taper(0.05, type='hann')
                tr.filter('lowpass',
                          zerophase=True,
                          freq=PAR.FREQ)

        elif PAR.FILTER == 'Highpass':
            for tr in traces:
                tr.detrend('demean')
                tr.detrend('linear')
                tr.taper(0.05, type='hann')
                tr.filter('highpass',
                          zerophase=True,
                          freq=PAR.FREQ)

        else:
            raise ParameterError()

        return traces

    def apply_mute(self, traces):
        if not PAR.MUTE:
            return traces

        if 'MuteEarlyArrivals' in PAR.MUTE:
            traces = signal.mute_early_arrivals(traces,
                PAR.MUTE_EARLY_ARRIVALS_SLOPE,  # (units: time/distance)
                PAR.MUTE_EARLY_ARRIVALS_CONST,  # (units: time)
                self.get_time_scheme(traces),
                self.get_source_coords(traces),
                self.get_receiver_coords(traces))

        if 'MuteLateArrivals' in PAR.MUTE:
            traces = signal.mute_late_arrivals(traces,
                PAR.MUTE_LATE_ARRIVALS_SLOPE,  # (units: time/distance)
                PAR.MUTE_LATE_ARRIVALS_CONST,  # (units: time)
                self.get_time_scheme(traces),
                self.get_source_coords(traces),
                self.get_receiver_coords(traces))

        if 'MuteShortOffsets' in PAR.MUTE:
            traces = signal.mute_short_offsets(traces,
                PAR.MUTE_SHORT_OFFSETS_DIST,
                self.get_source_coords(traces),
                self.get_receiver_coords(traces))

        if 'MuteLongOffsets' in PAR.MUTE:
            traces = signal.mute_long_offsets(traces,
                PAR.MUTE_LONG_OFFSETS_DIST,
                self.get_source_coords(traces),
                self.get_receiver_coords(traces))

        return traces

    def apply_normalize(self, traces):
        if not PAR.NORMALIZE:
            return traces

        if 'NormalizeEventsL1' in PAR.NORMALIZE:
            # normalize event by L1 norm of all traces
            w = 0.
            for tr in traces:
                w += np.linalg.norm(tr.data, ord=1)
            for tr in traces:
                tr.data /= w

        elif 'NormalizeEventsL2' in PAR.NORMALIZE:
            # normalize event by L2 norm of all traces
            w = 0.
            for tr in traces:
                w += np.linalg.norm(tr.data, ord=2)
            for tr in traces:
                tr.data /= w

        if 'NormalizeTracesL1' in PAR.NORMALIZE:
            # normalize each trace by its L1 norm
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=1)
                if w > 0:
                    tr.data /= w

        elif 'NormalizeTracesL2' in PAR.NORMALIZE:
            # normalize each trace by its L2 norm
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=2)
                if w > 0:
                    tr.data /= w

        return traces

    def apply_filter_backwards(self, traces):
        for tr in traces:
            tr.data = np.flip(tr.data)

        traces = self.apply_filter()

        for tr in traces:
            tr.data = np.flip(tr.data)

        return traces

    # Additional parameter checking
    def check_filter(self):
        """ Checks filter settings
        """
        assert getset(PAR.FILTER) < set([
            'Bandpass',
            'Lowpass',
            'Highpass'])

        if PAR.FILTER == 'Bandpass':
            if 'FREQMIN' not in PAR:
                raise ParameterError('FREQMIN')
            if 'FREQMAX' not in PAR:
                raise ParameterError('FREQMAX')
            assert 0 < PAR.FREQMIN
            assert PAR.FREQMIN < PAR.FREQMAX
            assert PAR.FREQMAX < np.inf

        elif PAR.FILTER == 'Lowpass':
            raise NotImplementedError
            if 'FREQ' not in PAR:
                raise ParameterError('FREQ')
            assert 0 < PAR.FREQ <= np.inf

        elif PAR.FILTER == 'Highpass':
            raise NotImplementedError
            if 'FREQ' not in PAR:
                raise ParameterError('FREQ')
            assert 0 <= PAR.FREQ < np.inf

    def check_mute(self):
        """ Checks mute settings
        """
        if not PAR.MUTE:
            return

        assert getset(PAR.MUTE) <= set([
            'MuteEarlyArrivals',
            'MuteLateArrivals',
            'MuteShortOffsets',
            'MuteLongOffsets'])

        if 'MuteEarlyArrivals' in PAR.MUTE:
            assert 'MUTE_EARLY_ARRIVALS_SLOPE' in PAR
            assert 'MUTE_EARLY_ARRIVALS_CONST' in PAR
            assert PAR.MUTE_EARLY_ARRIVALS_SLOPE >= 0.

        if 'MuteLateArrivals' in PAR.MUTE:
            assert 'MUTE_LATE_ARRIVALS_SLOPE' in PAR
            assert 'MUTE_LATE_ARRIVALS_CONST' in PAR
            assert PAR.MUTE_LATE_ARRIVALS_SLOPE >= 0.

        if 'MuteShortOffsets' in PAR.MUTE:
            assert 'MUTE_SHORT_OFFSETS_DIST' in PAR
            assert 0 < PAR.MUTE_SHORT_OFFSETS_DIST

        if 'MuteLongOffsets' in PAR.MUTE:
            assert 'MUTE_LONG_OFFSETS_DIST' in PAR
            assert 0 < PAR.MUTE_LONG_OFFSETS_DIST

        if 'MuteShortOffsets' not in PAR.MUTE:
            setattr(PAR, 'MUTE_SHORT_OFFSETS_DIST', 0.)

        if 'MuteLongOffsets' not in PAR.MUTE:
            setattr(PAR, 'MUTE_LONG_OFFSETS_DIST', 0.)

    def check_normalize(self):
        assert getset(PAR.NORMALIZE) < set([
            'NormalizeTracesL1',
            'NormalizeTracesL2',
            'NormalizeEventsL1',
            'NormalizeEventsL2'])

    # Utility functions
    def get_time_scheme(self, traces):
        """ FIXME: extract time scheme from trace headers rather than
            parameters file.
            Note from Alexis Bottero : it is actually better like this in
            my opinion because this allows for longer traces to be processed.
            Indeed, in su format only 2 bytes are dedicated to the number of
            samples which is supposed to be stored as an unsigned int. The
            maximum NT which can be stored in the header is then 32762 whereas
            there is no limit in principle.
        """
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
