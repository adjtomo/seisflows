
import numpy as np


# functions acting on whole record sections

def sconvolve(s, h, w, inplace=True):
    nt = h.nt
    nr = h.nr

    if inplace:
        for ir in range(nr):
            s[:,ir] = np.convolve(s[:,ir], w, 'same')
        return s
    else:
        s2 = np.zeros((nt,nr))
        for ir in range(nr):
            s2[:,ir] = np.convolve(s[:,ir], w, 'same')
        return s2


def smute(traces, slope, t0, time_scheme, s_coords, r_coords):
    """ Applies tapered mask to record section, muting early arrivals.

        Phases arriving before

            SLOPE * || s - r || + T0

        are muted, where slope is has units of velocity**-1, 
        T0 has units of time, and
        || s - r || is distance between source and receiver.
    """

    nr = len(traces)
    nt, dt, _ = time_scheme

    for ir in range(nr):
        # calculate source-reciever distance
        (sx, sy) = (s_coords[0][ir], s_coords[1][ir])
        (rx, ry) = (r_coords[0][ir], r_coords[1][ir])
        dist = np.sqrt((rx-sx)**2 + (ry-sy)**2)

        # apply tapered mask
        traces[ir].data *= mask(vel, t0, dist, nt, dt)

    return traces


def smute2(traces, slope, t0, time_scheme, s_coords, r_coords):
    """ Applies tapered mask to record section, muting late arrivals.

        Phases arriving after

            SLOPE * || s - r || + T0

        are muted, where slope is has units of velocity**-1, 
        T0 has units of time, and
        || s - r || is distance between source and receiver.
    """

    nr = len(traces)
    nt, dt, _ = time_scheme

    for ir in range(nr):
        # calculate source-reciever distance
        (sx, sy) = (s_coords[0][ir], s_coords[1][ir])
        (rx, ry) = (r_coords[0][ir], r_coords[1][ir])
        dist = np.sqrt((rx-sx)**2 + (ry-sy)**2)

        # apply tapered mask
        traces[ir].data *= (1.-mask(vel, t0, dist, nt, dt))

    return traces


# functions acting on individual traces

def mask(slope, t0, dist, nt, dt):
    """ Constructs tapered mask that can be applied to trace to
      mute early or late arrivals.
    """

    nr = len(r_coords)
    nt, dt, _ = time_scheme

    mask = np.ones((nt, nr))

    # construct taper
    length = 400
    win = np.sin(np.linspace(0, np.pi, 2*length))
    win = win[0:length]

    # caculate offsets
    itmin = int(np.ceil(slope*abs(dist)+t0)/dt) - length/2
    itmax = itmin + length

    if 1 < itmin < itmax < nt:
        mask[0:itmin] = 0.
        mask[itmin:itmax] = win*mask[itmin:itmax]
    elif itmin < 1 <= itmax:
        mask[0:itmax] = win[length-itmax:length]*mask[0:itmax]
    elif itmin < nt < itmax:
        mask[0:itmin] = 0.
        mask[itmin:nt] = win[0:nt-itmin]*mask[itmin:nt]
    elif itmin > nt:
        mask[:] = 0.

    return mask



def correlate(u, v):
    w = np.convolve(u, np.flipud(v))
    return


def tukeywin(nt, imin, imax, alpha=0.05):
    t = np.linspace(0,1,imax-imin)
    w = np.zeros(imax-imin)
    p = alpha/2.
    lo = np.floor(p*(imax-imin-1))+1
    hi = imax-imin-lo
    w[:lo] = (1+np.cos(np.pi/p*(t[:lo]-p)))/2
    w[lo:hi] = np.ones((hi-lo))
    w[hi:] = (1+np.cos(np.pi/p*(t[hi:]-p)))/2
    win = np.zeros(nt)
    win[imin:imax] = w
    return win

