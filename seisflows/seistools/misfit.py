
import numpy as np
from scipy.signal import hilbert


def Traveltime(wsyn, wobs, nt, dt):
    # cross correlation time
    cc = abs(np.convolve(wobs, np.flipud(wsyn)))
    cmax = 0
    misfit = 0.
    ioff = None
    for it in range(2*nt-1):
        if cc[it] > cmax:
            cmax = cc[it]
            ioff = it
            misfit = (ioff-nt+1)*dt
    if ioff is not None:
        misfit = (ioff-nt+1)*dt
    return misfit


def Amplitude(wsyn, wobs, nt, dt):
    # cross correlation amplitude
    cc = np.convolve(wobs, np.flipud(wsyn))
    cmax = 0
    ioff = 0
    for it in range(2*nt-1):
        if cc[it] > cmax:
            cmax = cc[it]
            ioff = it
    if ioff <= 0:
        wrsd = wsyn[ioff:] - wobs[:-ioff]
    else:
        wrsd = wsyn[:-ioff] - wobs[ioff:]
    return np.sqrt(np.sum(wrsd*wrsd*dt))


def Waveform(wsyn, wobs, nt, dt):
    # waveform rsderence
    wrsd = wsyn-wobs
    return np.sqrt(np.sum(wrsd*wrsd*dt))


def Envelope(wsyn, wobs, nt, dt, eps=0.05):
    # envelope rsderence
    esyn = abs(hilbert(wsyn))
    eobs = abs(hilbert(wobs))
    ersd = esyn-eobs
    return np.sqrt(np.sum(ersd*ersd*dt))


def InstantaneousPhase(wsyn, wobs, nt, dt, eps=0.05):
    # instantaneous phase 
    r = np.real(hilbert(wsyn))
    i = np.imag(hilbert(wsyn))
    phi_syn = np.arctan2(i,r)

    r = np.real(hilbert(wobs))
    i = np.imag(hilbert(wobs))
    phi_obs = np.arctan2(i,r)

    phi_rsd = phi_syn - phi_obs
    return np.sqrt(np.sum(phi_rsd*phi_rsd*dt))

