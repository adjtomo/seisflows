
import numpy as np
import obspy

from seisflows.tools import msg, unix
from seisflows.tools.code import exists, Struct
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

from seisflows.seistools import adjoint, misfit, readers, writers
from seisflows.seistools.signal import smute, smute2

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


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
        import solver
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


    def write_residuals(self, path, s, d):
        """ Computes residuals from observations and synthetics
        """
        nt, dt, _ = self.get_time_scheme(s)
        n, _ = self.get_network_size(s)

        filename = path +'/'+ 'residuals'
        if exists(filename):
            r = list(np.loadtxt(filename))
        else:
            r = []

        for i in range(n):
            r.append(self.misfit(s[i].data, d[i].data, nt, dt))

        np.savetxt(filename, r)


    def write_adjoint_traces(self, path, s, d, channel):
        """ Generates adjoint traces from observed and synthetic traces
        """
        nt, dt, _ = self.get_time_scheme(s)
        n, _ = self.get_network_size(s)

        for i in range(n):
            s[i].data = self.adjoint(s[i].data, d[i].data, nt, dt)

        self.writer(s, path, channel)


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

        elif PAR.MUTE in ['Linear', 'Early']:
            # mutes early arrivals
            return smute(traces, 
                PAR.MUTESLOPE, # (units: time/distance)
                PAR.MUTECONST, # (units: time)
                self.get_time_scheme(traces),
                self.get_source_coords(traces),
                self.get_receiver_coords(traces))

        elif PAR.MUTE in ['Late']:
            # mutes late arrivals
            return smute2(traces,
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
            # normalize each trace by its own power
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=1)
                if w > 0:
                    tr.data /= w
            return traces

        elif PAR.NORMALIZE == 'L2':
            # normalize each trace by its own power
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=2)
                if w > 0:
                    tr.data /= w
            return traces

        elif PAR.NORMALIZE == 'L2_squared':
            # normalize each trace by its own power
            for tr in traces:
                w = np.linalg.norm(tr.data, ord=2)
                if w > 0:
                    tr.data /= w**2.
            return traces


        elif PAR.NORMALIZE == 'L2_summed':
            raise NotImplementedError


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
            if 'FREQLO' not in PAR: raise ParameterError('FREQLO')
            if 'FREQHI' not in PAR: raise ParameterError('FREQHI')
            assert 0 < PAR.FREQLO
            assert PAR.FREQLO < PAR.FREQHI
            assert PAR.FREQHI < infinity

        elif PAR.FILTER == 'Lowpass':
            raise NotImplementedError
            if 'FREQ' not in PAR: raise ParameterError('FREQ')
            assert 0 < PAR.FREQHI <= infinity

        elif PAR.FILTER == 'Highpass':
            raise NotImplementedError
            if 'FREQ' not in PAR: raise ParameterError('FREQ')
            assert 0 <= PAR.FREQLO < infinity

        else:
            raise ParameterError()


    def check_mute(self):
        """ Checks mute settings
        """
        if not PAR.MUTE:
            pass

        elif PAR.MUTE in ['Early', 'Late', 'Linear']:
            assert 'MUTESLOPE' in PAR
            assert 'MUTECONST' in PAR
            assert PAR.MUTESLOPE >= 0.

        else:
            raise ParameterError()


    def check_normalize(self):
        if not PAR.FILTER:
            pass
        elif PAR.FILTER in ['L1', 'L2', 'L2_squared']:
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



