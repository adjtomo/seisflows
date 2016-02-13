
import numpy as np
from scipy.signal import hilbert

import misfit


### adjoint traces generators


def Traveltime(wsyn, wobs, nt, dt):
    # cross correlation traveltime
    # (Tromp et al. 2005, eq 45)
    wadj = np.zeros(nt)
    wadj[1:-1] = (wsyn[2:] - wsyn[0:-2])/(2.*dt)
    wadj *= 1./(sum(wadj*wadj)*dt)
    wadj *= misfit.wtime(wsyn,wobs,nt,dt)
    return wadj


def Amplitude(wsyn, wobs, nt, dt):
    # cross correlation amplitude
    wadj = 1./(sum(wsyn*wsyn)*dt) * wsyn
    wadj *= misfit.wampl(wsyn,wobs,nt,dt)
    return wadj


def Waveform(wsyn, wobs, nt, dt):
    # waveform rsderence
    # (Tromp et al. 2005, eq 9)
    wadj = wsyn - wobs
    return wadj


def Envelope(wsyn, wobs, nt, dt, eps=0.05):
    # envelope rsderence
    esyn = abs(hilbert(wsyn))
    eobs = abs(hilbert(wobs))
    etmp = (esyn - eobs)/(esyn + eps*esyn.max())
    wadj = etmp*wsyn - np.imag(hilbert(etmp*np.imag(hilbert(wsyn))))
    return wadj


def InstantaneousPhase(wsyn, wobs, nt, dt, eps=0.05):
    # instantaneous phase 
    r = np.real(hilbert(wsyn))
    i = np.imag(hilbert(wsyn))
    phi_syn = np.arctan2(i,r)

    r = np.real(hilbert(wobs))
    i = np.imag(hilbert(wobs))
    phi_obs = np.arctan2(i,r)

    phi_rsd = phi_syn - phi_obs
    esyn = abs(hilbert(wsyn))
    emax = max(esyn)

    wadj = phi_rsd*np.imag(hilbert(wsyn))/(esyn**2. + eps*emax) + \
           np.imag(hilbert(phi_rsd * wsyn)/(esyn**2. + eps*emax))

    return wadj


###

def precond1(d, s, h):
    s[1:-1,:] = (s[2:,:] - s[0:-2,:])/(2.*h.dt)
    s[1:-1,:] *= 1./(np.sum(s[1:-1,:]**2,axis=0)*h.dt)
    s[0,:] = 0.
    s[-1,:] = 0.
    return s


def precond2(d, s, h):
    return -s

