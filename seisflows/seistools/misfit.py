
import numpy as np
import scipy.signal


def traveltime(wsyn, wobs, nt, dt):
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


def amplitude(wsyn, wobs, nt, dt):
    # cross correlation amplitude
    cc = np.convolve(wobs, np.flipud(wsyn))
    cmax = 0
    ioff = 0
    for it in range(2*nt-1):
        if cc[it] > cmax:
            cmax = cc[it]
            ioff = it
    if ioff <= 0:
        wdiff = wsyn[ioff:] - wobs[:-ioff]
    else:
        wdiff = wsyn[:-ioff] - wobs[ioff:]
    return np.sqrt(np.sum(wdiff*wdiff*dt))


def waveform(wsyn, wobs, nt, dt):
    # waveform difference
    wdiff = wsyn-wobs
    return np.sqrt(np.sum(wdiff*wdiff*dt))


def envelope(wsyn, wobs, nt, dt, eps=0.05):
    # envelope difference
    esyn = abs(scipy.signal.hilbert(wsyn))
    eobs = abs(scipy.signal.hilbert(wobs))
    ediff = esyn-eobs
    return np.sqrt(np.sum(ediff*ediff*dt))


