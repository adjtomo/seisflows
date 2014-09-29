
import numpy as np

from seisflows.tools import unix
from seisflows.tools.codetools import Struct
from seisflows.tools.configure import getclass, ParameterObject
from seisflows.seistools import adjoint, misfit, sbandpass, smute

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')


class default(object):
  """ Data processing class
  """

  def __init__(self,reader=None,writer=None):
      """ Class constructor
      """

      # check user supplied parameters
      if 'MISFIT' not in PAR:
          raise Exception

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

      if 'MUTE' not in PAR:
	  setattr(PAR,'MUTE',False)

      if 'MUTESLOPE' not in PAR:
	  setattr(PAR,'MUTESLOPE',0.)

      if 'MUTECONST' not in PAR:
	  setattr(PAR,'MUTECONST',0.)

      if 'XCOMP' not in PAR:
	  setattr(PAR,'XCOMP',False)

      if 'YCOMP' not in PAR:
          setattr(PAR,'YCOMP',True)

      if 'ZCOMP' not in PAR:
          setattr(PAR,'ZCOMP',False)

      # IO routines
      self.reader = reader
      self.writer = writer


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
        system = getclass('system',PAR.SYSTEM)()
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
	# write adjoint traces needed for action of Jacobian
	d,h = self.load(prefix='traces/lcg')
	s,_ = self.load(prefix='traces/syn')

	d = self.apply(self.process_traces,[d],[h])
	s = self.apply(self.process_traces,[s],[h])

	s = self.apply(self.compute_adjoint,[s,d],[h,output_type])
	self.save(s,h)


      elif output_type == 4:
	# write adjoint traces needed for action of Hessian
	d,h = self.load(prefix='traces/obs')
	s,_ = self.load(prefix='traces/syn')

	d = self.apply(self.process_traces,[d],[h])
	s = self.apply(self.process_traces,[s],[h])

	s = self.apply(self.compute_adjoint,[s,d],[h,output_type])
	self.save(s,h)


      elif output_type == 5:
	# write adjoint traces needed for Hessian preconditioner
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
	  return s-s0

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

      # load data
      if PAR.XCOMP: (f.x,h.x) = self.reader(prefix=prefix,channel=1)
      if PAR.YCOMP: (f.y,h.y) = self.reader(prefix=prefix,channel=2)
      if PAR.ZCOMP: (f.z,h.z) = self.reader(prefix=prefix,channel=3)

      # check headers
      h = self.check_headers(h)

      return f,h


  def save(self,s,h,prefix='traces/adj/',suffix=''):
      """ Writes seismic data to disk
      """
      if PAR.XCOMP: self.writer(s.x,h,channel=1,prefix=prefix,suffix=suffix)
      if PAR.YCOMP: self.writer(s.y,h,channel=2,prefix=prefix,suffix=suffix)
      if PAR.ZCOMP: self.writer(s.z,h,channel=3,prefix=prefix,suffix=suffix)

      if not PAR.XCOMP: self.writer(np.zeros((h.nt,h.nr)),h,channel=1,
	  prefix=prefix,suffix=suffix)
      if not PAR.YCOMP: self.writer(np.zeros((h.nt,h.nr)),h,channel=2,
	  prefix=prefix,suffix=suffix)
      if not PAR.ZCOMP: self.writer(np.zeros((h.nt,h.nr)),h,channel=3,
	  prefix=prefix,suffix=suffix)


  def write_residuals(self,s,h):
      """ Writes residuals
      """
      # sum components
      sum = np.zeros((h.nr))
      if PAR.XCOMP: sum = sum + s.x
      if PAR.YCOMP: sum = sum + s.y
      if PAR.ZCOMP: sum = sum + s.z
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

      hdrs = headers.values()[1:]
      keys = ['dt','nt']
      for key in keys:
	for hdr in hdrs:
	  assert h[key] == hdr[key]

      return h

