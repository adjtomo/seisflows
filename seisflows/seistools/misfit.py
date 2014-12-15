
import numpy as np
import scipy.signal


def wtime(wsyn, wobs, nt, dt):
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


def wampl(wsyn, wobs, nt, dt):
    # cross correlation amplitude
    # FIXME: unitialized ioff, missing ioff == 0 case
    cc = np.convolve(wobs, np.flipud(wsyn))
    cmax = 0
    for it in range(2*nt-1):
        if cc[it] > cmax:
            cmax = cc[it]
            ioff = it
    try:
        if ioff < 0:
            wdiff = wsyn[ioff:]-wobs[:-ioff]
        if ioff > 0:
            wdiff = wsyn[:-ioff]-wobs[ioff:]
    except ArithmeticError:
        wdiff = wsyn-wobs
    return np.sqrt(np.sum(wdiff*wdiff*dt))


def wdiff(wsyn, wobs, nt, dt):
    # waveform difference
    wdiff = wsyn-wobs
    return np.sqrt(np.sum(wdiff*wdiff*dt))


def etime(wsyn, wobs, nt, dt):
    # envelope cross correlation
    pass


def ediff(wsyn, wobs, nt, dt, eps=0.05):
    # envelope difference
    esyn = abs(scipy.signal.hilbert(wsyn))
    eobs = abs(scipy.signal.hilbert(wobs))
    ediff = esyn-eobs
    return np.sqrt(np.sum(ediff*ediff*dt))


def cdiff(wsyn, wobs, nt, dt):
    cdiff = np.correlate(wobs, wsyn) - np.correlate(wobs, wobs)
    return np.sqrt(np.sum(cdiff*cdiff*dt))
