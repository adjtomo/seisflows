
import numpy as np
import obspy

from seisflows.tools import msg, unix
from seisflows.tools.code import exists, Struct, string_types
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
            setattr(PAR, 'MISFIT', 'Waveform')

        if 'CHANNELS' not in PAR:
            raise ParameterError(PAR, 'CHANNELS')

        if 'READER' not in PAR:
            raise ParameterError(PAR, 'READER')

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
        if PAR.READER not in dir(readers):
            print msg.ReaderError
            raise ParameterError()

        if PAR.WRITER not in dir(writers):
            print msg.WriterError
            raise ParameterError()

        if type(PAR.CHANNELS) not in string_types:
            print msg.ChannelError
            raise ParameterError()

        if PAR.READER != PAR.WRITER:
           print msg.DataFormatWarning % (PAR.READER, PAR.WRITER)


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
        for channel in self.channels:
            obs = self.reader(path+'/'+'traces/obs', channel)
            syn = self.reader(path+'/'+'traces/syn', channel)

            obs = self.process_traces(obs)
            syn = self.process_traces(syn)

            self.write_residuals(path, syn, obs)
            self.write_adjoint_traces(path+'/'+'traces/adj', syn, obs, channel)


    def process_traces(self, stream):
        """ Performs data processing operations on traces
        """
        print ' Warning bandpass filter disabled'

        nt, dt, _ = self.get_time_scheme(stream)
        n, _ = self.get_network_size(stream)
        df = dt**-1

        #for ir in range(n):
            #stream[ir].detrend()

            # filter data
            #if PAR.FREQLO and PAR.FREQHI:
            #    stream[ir].filter('bandpass', freqmin=PAR.FREQLO, freqmax=PAR.FREQHI)

        return stream


    def process_traces_adjoint(self, stream):
        """ Performs data processing operations on traces
        """
        print ' Warning bandpass filter disabled'

        nt, dt, _ = self.get_time_scheme(stream)
        n, _ = self.get_network_size(stream)
        df = dt**-1

        #for ir in range(n):
            #stream[ir].detrend()

            # filter data
            #if PAR.FREQLO and PAR.FREQHI:
            #    stream[ir].filter('bandpass', freqmin=PAR.FREQLO, freqmax=PAR.FREQHI)

        return stream



    def write_residuals(self, path, s, d):
        """ Computes residuals from observations and synthetics
        """
        nt, dt, _ = self.get_time_scheme(s)
        n, _ = self.get_network_size(s)

        filename = path +'/'+ 'residuals'
        if exists(filename):
            r = np.loadtxt(filename)
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

        # apply adjoint of filter
        #s = adjoint_filter(self.process_traces)(s)

        # normalize traces
        if PAR.NORMALIZE:
            for ir in range(n):
                w = np.linalg.norm(d[i], ord=2)
                if w > 0: 
                    s[:,ir] /= w

        # mute direct arrival
        if PAR.MUTE:
            vel = PAR.MUTESLOPE
            off = PAR.MUTECONST
            s = smute(s, h, vel, off, constant_spacing=False)

        self.writer(s, path, channel)


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

    def check_headers(self, data):
        """ Checks headers for consistency
        """
        nt, dt, _ = self.get_time_scheme(data)
        n, _ = self.get_network_size(data)

        if 'DT' in PAR:
            if dt != PAR.DT:
                print 'Warning: dt != PAR.DT'

        if 'NT' in PAR:
            if nt != PAR.NT:
                print 'Warning: nt != PAR.NT'

        if 'NREC' in PAR:
            if nr != PAR.NREC:
                print 'Warning: nr != PAR.NREC'

        return h


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


