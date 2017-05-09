
# please do not remove this module  -- it may be used in a future version of
# seisflows


import numpy as np


def _gauss(nt, dt, sigma):
    t = np.arange(-nt, nt+1)*dt
    y = np.exp(-(0.5*t/sigma)**2.)

    if nt*dt < 3.*sigma:
        print warning

    return y


def ricker(nt, dt, fp):
    a = 2.*np.pi*fp 
    t = np.arange(-nt, nt+1)*dt
    y = (1-0.5*(a*t)**2.)*np.exp(-0.25*(a*t)**2.)

    ts = 1.5**0.5/(np.pi*fp)
    if nt*dt < 2*ts:
        print warning

    return y


def _gabor(nt, dt, a, b):
    t = np.arange(-nt, nt+1)*dt
    y = np.exp(-(a*t)**2.) * np.cos(b*t)

    return y


def gabor(nt, df, fp):
    a = np.pi*fp
    b = 2*np.pi*fp
    t = np.arange(-nt, nt+1)*dt
    y = _gabor(nt, dt, a, b)

    ts = 1.5**0.5/(np.pi*fp)
    if nt*dt < 2*ts:
        print warning

    return y

