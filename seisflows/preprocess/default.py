
import numpy as np

from seisflows.tools import unix
from seisflows.tools.codetools import Struct
from seisflows.seistools import adjoint, misfit, sbandpass, smute
from seisflows.tools.configtools import ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class default(object):
    """ Data processing class
    """

    def check(self):
        """ Checks objects and parameters
        """

        # check objects
        if 'solver' not in OBJ:
            raise Excpetion

        if 'system' not in OBJ:
            raise Excpetion

        global solver
        import solver

        global system
        import system


        # check parameters
        if 'MISFIT' not in PAR:
            raise Exception

        if 'NORMALIZE' not in PAR:
            setattr(PAR,'NORMALIZE',True)


        # mute settings
        if 'MUTE' not in PAR:
            setattr(PAR,'MUTE',False)

        if 'MUTESLOPE' not in PAR:
            setattr(PAR,'MUTESLOPE',0.)

        if 'MUTECONST' not in PAR:
            setattr(PAR,'MUTECONST',0.)


        # filter settings
        if 'BANDPASS' not in PAR:
            setattr(PAR,'BANDPASS',False)

        if 'HIGHPASS' not in PAR:
            setattr(PAR,'HIGHPASS',False)

        if 'LOWPASS' not in PAR:
            setattr(PAR,'LOWPASS',False)

        if 'FREQLO' not in PAR:
            setattr(PAR,'FREQLO',0.)

        if 'FREQHI' not in PAR:
            setattr(PAR,'FREQHI',0.)


    def process_traces(self,s,h):
        """ Filters and mutes data
        """
        # filter data
        if PAR.BANDPASS:
            s = sbandpass(s,h,PAR.FREQLO,PAR.FREQHI)

        if PAR.HIGHPASS:
            s = shighpass(s,h,PAR.FREQLO)

        if PAR.HIGHPASS:
            s = slowpass(s,h,PAR.FREQHI)


        # mute direct arrival
        if PAR.MUTE == 1:
            vel = PAR.MUTESLOPE
            off = PAR.MUTECONST
            s = smute(s,h,vel,off,constant_spacing=False)

        elif PAR.MUTE == 2:
            vel = PAR.MUTESLOPE*(PAR.NREC+1)/(PAR.XMAX-PAR.XMIN)
            off = PAR.MUTECONST
            src = system.getnode()
            s = smute(s,h,vel,off,src,constant_spacing=True)

        return s


    def prepare_adjoint(self,path='.',output_type=2):
        """ Prepares solver for adjoint simulation by reading observations and
         synthetics, performing data processing, and writing various types of
         adjoint traces, depending on the keyword argument output_type.
        """
        unix.cd(path)

        if output_type == 1:
            # write residuals only
            d,h = self.load(prefix='traces/obs')
            s,_ = self.load(prefix='traces/syn')

            d = self.apply(self.process_traces,[d],[h])
            s = self.apply(self.process_traces,[s],[h])

            r = self.apply(self.compute_residuals,[s,d],[h],inplace=False)
            self.write_residuals(r,h)


        elif output_type == 2:
            # write adjoint traces needed for gradient evaluation
            d,h = self.load(prefix='traces/obs')
            s,_ = self.load(prefix='traces/syn')

            d = self.apply(self.process_traces,[d],[h])
            s = self.apply(self.process_traces,[s],[h])

            r = self.apply(self.compute_residuals,[s,d],[h],inplace=False)
            self.write_residuals(r,h)

            s = self.apply(self.compute_adjoint,[s,d],[h,output_type])
            self.save(s,h)


        elif output_type == 3:
            # write adjoint traces needed to get action of Jacobian
            d,h = self.load(prefix='traces/lcg')
            s,_ = self.load(prefix='traces/syn')

            d = self.apply(self.process_traces,[d],[h])
            s = self.apply(self.process_traces,[s],[h])

            s = self.apply(self.compute_adjoint,[s,d],[h,output_type])
            self.save(s,h)


        elif output_type == 4:
            # write adjoint traces needed to get action of Hessian
            d,h = self.load(prefix='traces/obs')
            s,_ = self.load(prefix='traces/syn')

            d = self.apply(self.process_traces,[d],[h])
            s = self.apply(self.process_traces,[s],[h])

            s = self.apply(self.compute_adjoint,[s,d],[h,output_type])
            self.save(s,h)


    def compute_residuals(self,s,d,h):
        """ Computes residuals from observations and synthetics
        """
        r = []
        for i in range(h.nr):
            r.append(self.call_misfit(s[:,i],d[:,i],h.nt,h.dt))
        return np.array(r)


    def compute_adjoint(self,s,d,h,output_type):
        """ Computes adjoint traces from observed and synthetic traces
        """
        if output_type == 1:
            pass

        elif output_type == 2:
            # gradient evaluation
            for i in range(h.nr):
                s[:,i] = self.call_adjoint(s[:,i],d[:,i],h.nt,h.dt)
            return s

        elif output_type == 3:
            # action of Jacobian
            return s-d

        elif output_type == 4:
            # action of Hessian
            pass

        elif output_type == 5:
            # Hessian preconditioner
            for i in range(h.nr):
                if PAR.HESSIAN == 1:
                    pass
                if PAR.HESSIAN == 2:
                    s[1:-1,i] = (s[2:,i]-2*s[1:-1,i]+s[0:-2,i]) / (h.dt**2)
                elif PAR.HESSIAN in [3,4]:
                    s[1:-1,i] = (s[2:,i]-s[0:-2,i]) / (2*h.dt)
                    s[:,i] = 1/(sum(f[:,i]*f[:,i])*h.dt) * f[:,i]
            return s

        # normalize
        if PAR.NORMALIZE==0:
            pass
        elif PAR.NORMALIZE==1:
            for ir in range(h.nr):
                s[:,ir] = s[:,ir]/np.norm(d[:,ir],ord=2)
        elif PAR.NORMALIZE==2:
            # normalize by trace
            for ir in range(h.nr):
                s[:,ir] = s[:,ir]/s[:,ir].max()
        elif PAR.NORMALIZE==3:
            # normalize by source
            s = s/s.max()

        return s


    ### misfit/adjoint wrappers

    def call_adjoint(self,wsyn,wobs,nt,dt):
        """ Wrapper for misfit functions
        """
        if PAR.MISFIT in ['wav','wdiff']:
            # waveform difference
            w = adjoint.wdiff(wsyn,wobs,nt,dt)
        elif PAR.MISFIT in ['tt','wtime']:
            # traveltime
            w = adjoint.wtime(wsyn,wobs,nt,dt)
        elif PAR.MISFIT in ['ampl','wampl']:
            # amplitude
            w = adjoint.wampl(wsyn,wobs,nt,dt)
        elif PAR.MISFIT in ['env','ediff']:
            # envelope
            w = adjoint.ediff(wsyn,wobs,nt,dt,eps=0.05)
        elif PAR.MISFIT in ['cdiff']:
            # cross correlation
            w = adjoint.cdiff(wsyn,wobs,nt,dt)
        else:
            w = wobs
        return w


    def call_misfit(self,wsyn,wobs,nt,dt):
        """ Wrapper for adjoint trace computations
        """
        if PAR.MISFIT in ['wav','wdiff']:
            # waveform difference
            e = misfit.wdiff(wsyn,wobs,nt,dt)
        elif PAR.MISFIT in ['tt','wtime']:
            # traveltime
            e = misfit.wtime(wsyn,wobs,nt,dt)
        elif PAR.MISFIT in ['ampl','wampl']:
            # amplitude
            e = misfit.wampl(wsyn,wobs,nt,dt)
        elif PAR.MISFIT in ['env','ediff']:
            # envelope
            e = misfit.ediff(wsyn,wobs,nt,dt,eps=0.05)
        elif PAR.MISFIT in ['cdiff']:
            # cross correlation
            e = misfit.cdiff(wsyn,wobs,nt,dt)
        else:
            e = 0.
        return float(e)



    ### input/output

    def load(self,prefix=''):
        """ Reads seismic data from disk
        """
        h = Struct()
        f = Struct()

        for channel in solver.channels:
            f[channel],h[channel] = solver.reader(prefix=prefix,channel=channel)

        # check headers
        h = self.check_headers(h)

        return f,h


    def save(self,s,h,prefix='traces/adj/',suffix=''):
        """ Writes seismic data to disk
        """
        for channel in solver.channels:
            solver.writer(s[channel],h,channel=channel,prefix=prefix,suffix=suffix)


    def write_residuals(self,s,h):
        """ Writes residuals
        """
        # sum components
        sum = np.zeros((h.nr))
        for channel in solver.channels:
            sum = sum + s[channel]
        np.savetxt('residuals',sum)



    ### utility functions

    def apply(self,func,arrays,args,inplace=True):
        """ Applies data processing operation to multi-component data
        """
        if inplace:
            if len(arrays) == 1:
                for key in arrays[0]:
                    arrays[0][key] = func(arrays[0][key],*args)
                return arrays[0]

            if len(arrays) == 2:
                for key in arrays[0]:
                    arrays[0][key] = func(arrays[0][key],arrays[1][key],*args)
                return arrays[0]

        else:
            new = Struct()
            if len(arrays) == 1:
                for key in arrays[0]:
                    new[key] = func(arrays[0][key],*args)
                return new

            if len(arrays) == 2:
                for key in arrays[0]:
                    new[key] = func(arrays[0][key],arrays[1][key],*args)
                return new


    def check_headers(self,headers):
        """ Checks headers for consistency
        """
        h = headers.values()[0]

        if h.dt != PAR.DT:
            h.dt = PAR.DT
        if h.nt != PAR.NT:
            print 'Warning: h.nt != PAR.NT'
        if h.nr != PAR.NREC:
            print 'Warning: h.nr != PAR.NREC'

        #hdrs = headers.values()[1:]
        #keys = ['dt','nt']
        #for key in keys:
        #  for hdr in hdrs:
        #    assert h[key] == hdr[key]

        return h


