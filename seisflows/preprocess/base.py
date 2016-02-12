
import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import Struct, string_types
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

from seisflows.seistools import adjoint, misfit, sbandpass, smute, readers, writers

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class base(object):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'MISFIT' not in PAR:
            setattr(PAR, 'MISFIT', 'waveform')

        if 'CHANNELS' not in PAR:
            raise ParameterError(PAR, 'CHANNELS')

        if 'READER' not in PAR:
            raise ParameterError(PAR, 'CHANNELS')

        if 'WRITER' not in PAR:
            setattr(PAR, 'WRITER', PAR.READER)

        if 'NORMALIZE' not in PAR:
            setattr(PAR, 'NORMALIZE', True)

        # mute settings
        if 'MUTE' not in PAR:
            setattr(PAR, 'MUTE', False)

        if 'MUTESLOPE' not in PAR:
            setattr(PAR, 'MUTESLOPE', 0.)

        if 'MUTECONST' not in PAR:
            setattr(PAR, 'MUTECONST', 0.)

        # filter settings
        if 'BANDPASS' not in PAR:
            setattr(PAR, 'BANDPASS', False)

        if 'FREQLO' not in PAR:
            setattr(PAR, 'FREQLO', 0.)

        if 'FREQHI' not in PAR:
            setattr(PAR, 'FREQHI', 0.)

        # assertions
        if PAR.READER.lower() not in dir(readers):
            #print msg.ReaderError
            raise ParameterError()

        if PAR.WRITER.lower() not in dir(writers):
            #print msg.WriterError
            raise ParameterError()

        if type(PAR.CHANNELS) not in string_types:
            #print msg.ChannelError
            raise ParameterError()


    def setup(self):
        """ Sets up data preprocessing machinery
        """
        # define misfit function
        self.misfit = getattr(misfit, PAR.MISFIT.lower())

        # define adjoint trace generator
        self.adjoint = getattr(adjoint, PAR.MISFIT.lower())

        # define seismic data reader
        self.reader = getattr(readers, PAR.READER)

        # define seismic data writer
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

        r = self.multichannel(self.write_residuals, [s, d], [h], inplace=False)
        s = self.multichannel(self.generate_adjoint_traces, [s, d], [h])
        self.save(s, h, prefix='traces/adj/')


    def process_traces(self, s, h):
        """ Performs data processing operations on traces
        """
        # filter data
        if PAR.BANDPASS:
            if PAR.FREQLO and PAR.FREQHI:
                s = sbandpass(s, h, PAR.FREQLO, PAR.FREQHI)

            elif PAR.FREQHI:
                s = shighpass(s, h, PAR.FREQLO)

            elif PAR.FREQHI:
                s = slowpass(s, h, PAR.FREQHI)

            else:
                raise ParameterError(PAR, 'BANDPASS')

        # mute direct arrival
        if PAR.MUTE:
            vel = PAR.MUTESLOPE
            off = PAR.MUTECONST
            s = smute(s, h, vel, off, constant_spacing=False)

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

        for i in range(h.nr):
            s[:,i] = self.adjoint(s[:,i], d[:,i], h.nt, h.dt)

        # bandpass once more
        if PAR.BANDPASS:
            if PAR.FREQLO and PAR.FREQHI:
                s = sbandpass(s, h, PAR.FREQLO, PAR.FREQHI, 'reverse')

            elif PAR.FREQHI:
                s = shighpass(s, h, PAR.FREQLO, 'reverse')

            elif PAR.FREQHI:
                s = slowpass(s, h, PAR.FREQHI, 'reverse')

            else:
                raise ParameterError(PAR, 'BANDPASS')

        # normalize traces
        if PAR.NORMALIZE:
            for ir in range(h.nr):
                w = np.linalg.norm(d[:,ir], ord=2)
                if w > 0: 
                    s[:,ir] /= w

        return s


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


