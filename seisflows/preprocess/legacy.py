
import numpy as np
import scipy.signal as _signal

from seisflows.tools import unix
from seisflows.tools.code import Struct
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError
from seisflows.tools.math import infinity

from seisflows.seistools import adjoint, misfit, readers, writers
from seisflows.seistools.signal import sbandpass, smute, smutelow

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class legacy(object):
    """ Legacy data preprocessing class

      Legacy methods differ from base methods mainly in that the former depend
      depend on scipy and the latter depend on obspy I/O and signal processing
      functions.
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'MISFIT' not in PAR:
            setattr(PAR, 'MISFIT', None)

        if 'BACKPROJECT' not in PAR:
            setattr(PAR, 'BACKPROJECT', None)

        if 'CHANNELS' not in PAR:
            raise ParameterError(PAR, 'CHANNELS')

        if 'READER' not in PAR:
            raise ParameterError(PAR, 'READER')

        if 'WRITER' not in PAR:
            setattr(PAR, 'WRITER', PAR.READER)

        if 'NORMALIZE' not in PAR:
            setattr(PAR, 'NORMALIZE', 'L2')

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
        # define misfit function and adjoint source generator
        if PAR.MISFIT:
            self.misfit = getattr(misfit, PAR.MISFIT)
            self.adjoint = getattr(adjoint, PAR.MISFIT)
        elif PAR.BACKPROJECT:
            self.adjoint = getattr(adjoint, PAR.BACKPROJECT)

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
        unix.cd(path)

        d, h = self.load(prefix='traces/obs/')
        s, _ = self.load(prefix='traces/syn/')

        d = self.multichannel(self.process_traces, [d], [h])
        s = self.multichannel(self.process_traces, [s], [h])

        s = self.multichannel(self.generate_adjoint_traces, [s, d], [h])
        self.save(s, h, prefix='traces/adj/')

        if PAR.MISFIT:
            self.multichannel(self.write_residuals, [s, d], [h], inplace=False)


    def process_traces(self, s, h):
        """ Performs data processing operations on traces
        """
        s = self.apply_filter(s, h)
        s = self.apply_normalize(s, h)
        s = self.apply_mute(s, h)

        return s


    def write_residuals(self, s, d, h):
        """ Computes residuals from observations and synthetics
        """
        r = []
        for i in range(h.nr):
            r.append(self.misfit(s[:,i], d[:,i], h.nt, h.dt))

        # write residuals
        np.savetxt('residuals', r)

        return np.array(r)


    def generate_adjoint_traces(self, s, d, h):
        """ Generates adjoint traces from observed and synthetic traces
        """
        # generate adjoint traces
        for i in range(h.nr):
            s[:,i] = self.adjoint(s[:,i], d[:,i], h.nt, h.dt)

        return s


    ### signal processing

    def apply_filter(self, s, h):
        if not PAR.FILTER:
            return s

        elif PAR.FILTER == 'Bandpass':
            s = _signal.detrend(s)
            return sbandpass(s, h, PAR.FREQLO, PAR.FREQHI)

        elif PAR.FILTER == 'Lowpass':
            raise NotImplementedError

        elif PAR.FILTER == 'Highpass':
            raise NotImplementedError

        elif PAR.FILTER == 'Butterworth':
            raise NotImplementedError

        else:
            raise ParameterError()


    def apply_mute(self, s, h):
        if not PAR.MUTE:
            return s

        elif PAR.MUTE == 'Simple':
            vel = PAR.MUTESLOPE
            off = PAR.MUTECONST

            # mute early arrivals
            return smute(s, h, vel, off, 0., constant_spacing=False)

        else:
            raise ParameterError()


    def apply_normalize(self, s, h):
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
            # normalize traces by their combined power
            w = np.linalg.norm(s, ord=2)
            if w > 0:
                s /= w
            return s

        else:
            raise ParameterError()


    ### input/output

    def load(self, prefix=None, suffix=None):
        """ Reads seismic data from disk
        """
        kwargs = dict()
        if prefix != None:
            kwargs['prefix'] = prefix
        if suffix != None:
            kwargs['suffix'] = suffix

        f = Struct()
        h = Struct()
        for channel in self.channels:
            f[channel], h[channel] = self.reader(channel=channel, **kwargs)

        # check headers
        h = self.check_headers(h)

        return f, h


    def save(self, s, h, prefix='traces/adj/', suffix=None):
        """ Writes seismic data to disk
        """
        kwargs = dict()
        if prefix != None:
            kwargs['prefix'] = prefix
        if suffix != None:
            kwargs['suffix'] = suffix

        for channel in self.channels:
            self.writer(s[channel], h, channel=channel, **kwargs)


    ### utility functions

    def multichannel(self, func, arrays, input, inplace=True):
        """ Applies function to multi-component data
        """
        if inplace:
            output = Struct(arrays[0])
        else:
            output = Struct()

        for channel in self.channels:
            args = [array[channel] for array in arrays] + input
            output[channel] = func(*args)

        return output


    def check_headers(self, headers):
        """ Checks headers for consistency
        """
        h = headers.values()[0]

        if 'DT' in PAR:
            if h.dt != PAR.DT:
                h.dt = PAR.DT

        if 'NT' in PAR:
            if h.nt != PAR.NT:
                print 'Warning: h.nt != PAR.NT'

        if 'NREC' in PAR:
            if h.nr != PAR.NREC:
                print 'Warning: h.nr != PAR.NREC'

        return h


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


    def get_size(self, channels):
        for channel in channels:
            _, h = self.reader(channel=channel, **kwargs)    
            print channel, h.nr, h.nt
