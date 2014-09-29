
import numpy as np

from seisflows.tools import signal


def sbandpass(s,h,freqlo,freqhi):
  nr = h.nr
  dt = h.dt
  fs = 1/dt

  for ir in range(0,nr):
    s[:,ir] = signal.bandpass(s[:,ir],freqlo,freqhi,fs)
  return s


def slowpass(s,h,freq):
   raise NotImplementedError


def shighpass(s,h,freq):
   raise NotImplementedError


def smute(s,h,vel,toff,xoff=0,constant_spacing=False):
  nt = h.nt
  dt = h.dt
  nr = h.nr

  # construct tapered window
  length = 400
  win = np.sin(np.linspace(0,1,2*length))
  win = win[0:length]

  for ir in range(0,nr):

    # calculate slope
    if vel!=0:
      slope = 1./vel
    elif vel==0:
      slope = 0

    # calculate offsets
    if constant_spacing:
      itoff = toff/dt
      ixoff = (ir-xoff)/dt
    else:
      itoff = toff/dt
      ixoff = (h.rx[ir]-h.sx[0]-xoff)/dt

    itmin = int(np.ceil(slope*abs(ixoff)+itoff)) - length/2
    itmax = itmin + length

    # apply window
    if 1 < itmin and itmax < nt:
      s[0:itmin,ir] = 0.
      s[itmin:itmax,ir] = win*s[itmin:itmax,ir]
    elif itmin < 1 and itmax >= 1:
      s[1:itmax,ir] = win[length-itmax+1:length]*s[1:itmax,ir]
    elif itmin < nt and itmax > nt:
      s[0:itmin,ir] = 0.
      s[itmin:nt,ir] = win[1:nt-itmin+1]*s[itmin:nt,ir]
    elif itmin > nt:
      s[:,ir] = 0.

  return s


def swindow(s,h,tmin,tmax,alpha=0.05,units='samples'):
    nt = h.nt
    dt = h.dt
    nr = h.nr
    t0 = h.t0

    if units == 'time':
      itmin = int((tmin-t0)/dt)
      itmax = int((tmax-t0)/dt)

    elif units == 'samples':
      itmin = int(tmin)
      itmax = int(tmax)

    win = signal.tukeywin(nt,itmin,itmax,alpha)

    # apply window
    for ir in range(0,nr):
      s[:,ir] = win*s[:,ir]

    return s


