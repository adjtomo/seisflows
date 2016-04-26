
import numpy as np
import obspy

from seisflows.tools import msg, unix
from seisflows.tools.code import exists, Struct
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

from seisflows.seistools import adjoint, misfit, readers, writers
from seisflows.seistools.signal import smute

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class base(object):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'MISFIT' not in PAR:
            setattr(PAR, 'MISFIT', 'Waveform')

        if 'CHANNELS' not in PAR:
            raise ParameterError(PAR, 'CHANNELS')

        if 'READER' not in PAR:
            raise ParameterError(PAR, 'READER')

        if 'WRITER' not in PAR:
            setattr(PAR, 'WRITER', PAR.READER)

        if 'NORMALIZE' not in PAR:
            setattr(PAR, 'NORMALIZE', True)

        if 'MUTE' not in PAR:
            setattr(PAR, 'MUTE', None)

        if 'FILTER' not in PAR:
            setattr(PAR, 'FILTER', None)

        # assertions
        if PAR.READER not in dir(readers):
            print msg.ReaderError
            raise ParameterError()

        if PAR.WRITER not in dir(writers):
            print msg.WriterError
            raise ParameterError()

        self.check_filter()
        self.check_mute()
        self.check_normalize()


    def setup(self):
        """ Sets up data preprocessing machinery
        """
        # define misfit function and adjoint trace generator
        self.misfit = getattr(misfit, PAR.MISFIT)
        self.adjoint = getattr(adjoint, PAR.MISFIT)

        # define seismic data reader and writer
        self.reader = getattr(readers, PAR.READER)
        self.writer = getattr(writers, PAR.WRITER)

        # prepare channels list
        self.channels = []
        for char in PAR.CHANNELS:
            self.channels += [char]


    def prepare_eval_grad(self, path='.'):
        """ Prepares solver for gradient evaluation by writing residuals and
          adjoint traces
        """
        import solver
        for channel in solver.channels:
            obs = self.reader(path+'/'+'traces/obs', channel)
            syn = self.reader(path+'/'+'traces/syn', channel)

            # process observations
            obs = self.apply_filter(obs)
            obs = self.apply_mute(obs)
            obs = self.apply_normalize(obs)

            # process synthetics
            syn = self.apply_filter(syn)
            syn = self.apply_mute(syn)
            syn = self.apply_normalize(syn)

            self.write_residuals(path, syn, obs)
            self.write_adjoint_traces(path+'/'+'traces/adj', syn, obs, channel)


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

        # generate adjoint traces
        for i in range(n):
            s[i].data = self.adjoint(s[i].data, d[i].data, nt, dt)

        self.writer(s, path, channel)


    ### signal processing

    def apply_filter(self, s):
        if not PAR.FILTER:
            return s

        elif PAR.FILTER == 'Simple':
            s = _signal.detrend(s)
            for trace in s:
                trace.filter('bandpass', freqmin=PAR.FREQLO, freqmax=PAR.FREQHI)

        elif PAR.FILTER == 'Butterworth':
            raise NotImplementedError

        else:
            raise ParameterError()

        # workaround obspy dtype conversion
        trace.data = trace.data.astype(np.float32)


    def apply_mute(self, s):
        if not PAR.MUTE:
            return s

        elif PAR.MUTE == 'Simple':
            vel = PAR.MUTESLOPE
            off = PAR.MUTECONST

            # mute early arrivals
            return smute(s, h, vel, off, 0., constant_spacing=False)

        else:
            raise ParameterError()


    def apply_normalize(self, s):
        if not PAR.NORMALIZE:
            return s

        elif PAR.NORMALIZE == 'L2':
            # normalize each trace by its own power
            for ir in range(h.nr):
                w = np.linalg.norm(s[:,ir], ord=2)
                if w > 0:
                    s[:,ir] /= w
            return s

        elif PAR.NORMALIZE == 'L2_all':
            # normalize all traces by their combined power
            w = np.linalg.norm(d, ord=2)
            if w > 0:
                s /= w
            return s


    def apply_filter_backwards(self, s):
        for trace in s:
            trace.data = np.flip(trace.data)

        s = self.apply_filter()

        for trace in s:
            trace.data = np.flip(trace.data)

        return s



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
            if 'FREQLO' not in PAR: raise ParameterError('FREQLO')
            if 'FREQHI' not in PAR: raise ParameterError('FREQHI')
            assert 0 == PAR.FREQLO
            assert PAR.FREQHI <= infinity

        elif PAR.FILTER == 'Highpass':
            raise NotImplementedError
            if 'FREQLO' not in PAR: raise ParameterError('FREQLO')
            if 'FREQHI' not in PAR: raise ParameterError('FREQHI')
            assert 0 <= PAR.FREQLO
            assert PAR.FREQHI == infinity

        elif PAR.FILTER == 'Butterworth':
            raise NotImplementedError
            if 'CORNERS' not in PAR: raise ParameterError
            if 'NPASS' not in PAR: PAR.NPASS = 2
            assert len(PAR.CORNERS) == 4
            assert 0. <= PAR.CORNERS[1]
            assert PAR.CORNERS[0] < PAR.CORNERS[1]
            assert PAR.CORNERS[1] < PAR.CORNERS[2]
            assert PAR.CORNERS[2] < PAR.CORNERS[3]
            assert PAR.CORNERS[3] <= infinity

        else:
            raise ParameterError()


    def check_mute(self):
        """ Checks mute settings
        """
        if not PAR.MUTE:
            pass

        elif PAR.MUTE == 'Simple':
            if 'MUTESLOPE' not in PAR: PAR.MUTESLOPE = infinity
            if 'MUTECONST' not in PAR: PAR.MUTECONST = 0.
            assert PAR.MUTESLOPE > 0.
            assert PAR.MUTECONST >= 0.

        else:
            raise ParameterError()


    def check_normalize(self):
        pass



    ### utility functions

    def write_zero_traces(self, path, channel):
        from obspy.core.stream import Stream
        from obspy.core.trace import Trace

        # construct seismic data and headers
        t = Trace(data=np.zeros(PAR.NT, dtype='float32'))
        t.stats.delta = PAR.DT
        s = Stream(t)*PAR.NREC

        # write to disk
        self.writer(s, path, channel)


    def get_time_scheme(self, stream):
        dt = PAR.DT
        nt = PAR.NT
        t0 = 0.
        return nt, dt, t0


    def get_network_size(self, stream):
        nrec = len(stream)
        nsrc = 1
        return nrec, nsrc


    def get_receiver_coords(self, stream):
        xr = []
        yr = []
        zr = []
        return xr, yr, zr


    def get_source_coords(self, stream):
        xr = []
        yr = []
        zr = []
        return xs, ys, zs


