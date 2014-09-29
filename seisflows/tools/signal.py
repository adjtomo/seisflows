
import numpy as _np
import scipy.signal as _signal


def correlate(u,v):
    w = _np.convolve(u,_np.flipud(v))
    return


def bandpass(w,freqlo,freqhi,fs,npass=2):
    wn = [2*freqlo/fs,2*freqhi/fs]
    (b,a) = _signal.butter(npass,wn)
    w = _signal.filtfilt(b,a,w)
    return w


def highpass(w,freq,df,npass=2):
    raise Exception('Not yet implemented.')


def lowpass(w,freq,df,npass=2):
    raise Exception('Not yet implemented.')


def window(nt,type='sine',**kwargs):
    if type=='sine':
      return _np.sin(_np.linspace(0,1,2*length))
    elif type=='tukey':
      return tukeywin


def tukeywin(nt,imin,imax,alpha=0.05):
    t = _np.linspace(0,1,imax-imin)
    w = _np.zeros(imax-imin)
    p = alpha/2.
    lo = _np.floor(p*(imax-imin-1))+1
    hi = imax-imin-lo
    w[:lo] = (1+_np.cos(_np.pi/p*(t[:lo]-p)))/2
    w[lo:hi] = _np.ones((hi-lo))
    w[hi:] = (1+_np.cos(_np.pi/p*(t[hi:]-p)))/2
    win = _np.zeros(nt)
    win[imin:imax] = w
    return win

