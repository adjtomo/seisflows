
import numpy as _np
from scipy.signal import hilbert as _analytic


def Traveltime(syn, obs, nt, dt):
    # cross correlation time
    cc = abs(_np.convolve(obs, _np.flipud(syn)))
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


def Amplitude(syn, obs, nt, dt):
    # cross correlation amplitude
    cc = _np.convolve(obs, _np.flipud(syn))
    cmax = 0
    ioff = 0
    for it in range(2*nt-1):
        if cc[it] > cmax:
            cmax = cc[it]
            ioff = it
    if ioff <= 0:
        wrsd = syn[ioff:] - obs[:-ioff]
    else:
        wrsd = syn[:-ioff] - obs[ioff:]
    return _np.sqrt(_np.sum(wrsd*wrsd*dt))


def Waveform(syn, obs, nt, dt):
    # waveform difference
    wrsd = syn-obs
    return _np.sqrt(_np.sum(wrsd*wrsd*dt))


def Envelope(syn, obs, nt, dt, eps=0.05):
    # envelope difference
    # (Yuan et al 2015, eq 9)
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    ersd = esyn-eobs
    return _np.sqrt(_np.sum(ersd*ersd*dt))


def InstantaneousPhase(syn, obs, nt, dt, eps=0.05):
    # instantaneous phase 
    # from Bozdag et al. 2011

    r = _np.real(_analytic(syn))
    i = _np.imag(_analytic(syn))
    phi_syn = _np.arctan2(i,r)

    r = _np.real(_analytic(obs))
    i = _np.imag(_analytic(obs))
    phi_obs = _np.arctan2(i,r)

    phi_rsd = phi_syn - phi_obs
    return _np.sqrt(_np.sum(phi_rsd*phi_rsd*dt))


def Envelope2(syn, obs, nt, dt, eps=0.):
    # envelope amplitude ratio
    # (Yuan et al 2015, eq B-1)
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    raise NotImplementedError


def Envelope3(syn, obs, nt, dt, eps=0.):
    # envelope cross-correlation lag
    # (Yuan et al 2015, eqs B-4)
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    return Traveltime(esyn, eobs, nt, dt)


def AnalyticSignal(syn, obs, nt, dt, eps=0.):
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))

    esyn1 = esyn + eps*max(esyn)
    eobs1 = eobs + eps*max(eobs)

    diff = syn/esyn1 - obs/eobs1

    return _np.sqrt(_np.sum(diff*diff*dt))

