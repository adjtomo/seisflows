#
# This is Seisflows
#
# See LICENCE file
#
# Functions used by the PREPROCESS class and specified by the MISFIT parameter
###############################################################################

# Import Numpy and utilities from Scipy
import numpy as np
from scipy.signal import hilbert as _analytic


def Waveform(syn, obs, nt, dt):
    """ Waveform difference
    """
    wrsd = syn-obs
    return np.sqrt(np.sum(wrsd*wrsd*dt))


def Envelope(syn, obs, nt, dt, eps=0.05):
    """ Envelope difference (Yuan et al 2015, eq 9)
    """
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    ersd = esyn-eobs
    return np.sqrt(np.sum(ersd*ersd*dt))


def InstantaneousPhase(syn, obs, nt, dt, eps=0.05):
    """ Instantaneous phase from Bozdag et al. 2011
    """
    r = np.real(_analytic(syn))
    i = np.imag(_analytic(syn))
    phi_syn = np.arctan2(i, r)

    r = np.real(_analytic(obs))
    i = np.imag(_analytic(obs))
    phi_obs = np.arctan2(i, r)

    phi_rsd = phi_syn - phi_obs
    return np.sqrt(np.sum(phi_rsd*phi_rsd*dt))


def Traveltime(syn, obs, nt, dt):
    """ Compute cross correlation traveltime between two traces suposing that
        they contain only one arrival
    """
    cc = abs(np.convolve(obs, np.flipud(syn)))
    return (np.argmax(cc)-nt+1)*dt


def TraveltimeInexact(syn, obs, nt, dt):
    """ Much faster (but possibly inaccurate) version of Traveltime function
    """
    it = np.argmax(syn)
    jt = np.argmax(obs)
    return (jt-it)*dt


def Amplitude(syn, obs, nt, dt):
    """ Cross correlation amplitude
    """
    ioff = (np.argmax(cc)-nt+1)*dt
    if ioff <= 0:
        wrsd = syn[ioff:] - obs[:-ioff]
    else:
        wrsd = syn[:-ioff] - obs[ioff:]
    return np.sqrt(np.sum(wrsd*wrsd*dt))


def Envelope2(syn, obs, nt, dt, eps=0.):
    """ Envelope amplitude ratio (Yuan et al 2015, eq B-1)
    """
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    raise NotImplementedError


def Envelope3(syn, obs, nt, dt, eps=0.):
    """ Envelope cross-correlation lag (Yuan et al 2015, eqs B-4)
    """
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))
    return Traveltime(esyn, eobs, nt, dt)


def InstantaneousPhase2(syn, obs, nt, dt, eps=0.):
    esyn = abs(_analytic(syn))
    eobs = abs(_analytic(obs))

    esyn1 = esyn + eps*max(esyn)
    eobs1 = eobs + eps*max(eobs)

    diff = syn/esyn1 - obs/eobs1

    return np.sqrt(np.sum(diff*diff*dt))


def Displacement(syn, obs, nt, dt):
    return Exception('This function can only used for migration.')

def Velocity(syn, obs, nt, dt):
    return Exception('This function can only used for migration.')

def Acceleration(syn, obs, nt, dt):
    return Exception('This function can only used for migration.')
