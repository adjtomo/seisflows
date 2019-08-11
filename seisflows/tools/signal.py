#
# This is Seisflows
#
# See LICENCE file
#
###############################################################################

# Import Numpy
import numpy as np


# Functions acting on whole record sections

def sconvolve(s, h, w, inplace=True):
    nt = h.nt
    nr = h.nr

    if inplace:
        for ir in range(nr):
            s[:, ir] = np.convolve(s[:, ir], w, 'same')
        return s
    else:
        s2 = np.zeros((nt, nr))
        for ir in range(nr):
            s2[:, ir] = np.convolve(s[:, ir], w, 'same')
        return s2


def mute_early_arrivals(traces, slope, const, time_scheme, s_coords, r_coords):
    """ Applies tapered mask to record section, muting early arrivals

        Signals arriving before

            SLOPE * || s - r || + CONST

        are muted, where slope is has units of velocity**-1,
        CONST has units of time, and
        || s - r || is distance between source and receiver.
    """

    nr = len(traces)
    nt, dt, _ = time_scheme

    for ir in range(nr):
        # calculate source-reciever distance
        (sx, sy) = (s_coords[0][ir], s_coords[1][ir])
        (rx, ry) = (r_coords[0][ir], r_coords[1][ir])
        offset = np.sqrt((rx-sx)**2 + (ry-sy)**2)

        # apply tapered mask
        traces[ir].data *= mask(slope, const, offset, (nt, dt, 0.))

    return traces


def mute_late_arrivals(traces, slope, const, time_scheme, s_coords, r_coords):
    """ Applies tapered mask to record section, muting late arrivals

        Signals arriving after

            SLOPE * || s - r || + CONST

        are muted, where SLOPE is has units of velocity**-1,
        CONST has units of time, and
        || s - r || is distance between source and receiver.
    """

    nr = len(traces)
    nt, dt, _ = time_scheme

    for ir in range(nr):
        # calculate source-reciever distance
        (sx, sy) = (s_coords[0][ir], s_coords[1][ir])
        (rx, ry) = (r_coords[0][ir], r_coords[1][ir])
        offset = np.sqrt((rx-sx)**2 + (ry-sy)**2)

        # apply tapered mask
        traces[ir].data *= (1.-mask(slope, const, offset, (nt, dt, 0.)))

    return traces


def mute_short_offsets(traces, dist, s_coords, r_coords):
    """ Mutes traces having

            || s - r || < DIST

        where || s - r || is the offset between source and receiver and
        DIST is a user-supplied cutoff
    """
    nr = len(traces)

    for ir in range(nr):
        # calculate source-reciever distance
        (sx, sy) = (s_coords[0][ir], s_coords[1][ir])
        (rx, ry) = (r_coords[0][ir], r_coords[1][ir])
        offset = np.sqrt((rx-sx)**2 + (ry-sy)**2)

        if offset < dist:
            traces[ir].data[:] = 0.

    return traces


def mute_long_offsets(traces, dist, s_coords, r_coords):
    """ Mutes traces having

            || s - r || > DIST

        where || s - r || is the offset between source and receiver and
        DIST is a user-supplied cutoff
    """
    nr = len(traces)

    for ir in range(nr):
        # calculate source-reciever distance
        (sx, sy) = (s_coords[0][ir], s_coords[1][ir])
        (rx, ry) = (r_coords[0][ir], r_coords[1][ir])
        offset = np.sqrt((rx-sx)**2 + (ry-sy)**2)

        if offset > dist:
            traces[ir].data[:] = 0.

    return traces


# Functions acting on individual traces

def mask(slope, const, offset, time_scheme, length=400):
    """ Constructs tapered mask that can be applied to trace to
      mute early or late arrivals.
    """

    nt, dt, _ = time_scheme

    mask = np.ones(nt)

    # construct taper
    win = np.sin(np.linspace(0, np.pi, 2*length))
    win = win[0:length]

    # caculate offsets
    itmin = int(np.ceil((slope*abs(offset)+const)/dt)) - length/2
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
    return w


def tukeywin(nt, imin, imax, alpha=0.05):
    t = np.linspace(0, 1, imax-imin)
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
