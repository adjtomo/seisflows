
import numpy as np


def ricker(dt, f0, width=100.):

    tb = 0.88521*dt
    nt = np.ceil(2*width*tb/dt)

    t = np.linspace(-width*tb, width*tb)
    a = (np.pi*f0)**2
    y = (1. - 2.*a*t**2) * np.exp(-a*t**2)

    return y


def gabor(nt, dt, f0):
    x = np.arange(-nt, nt)*dt


