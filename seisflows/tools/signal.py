"""
Signal processing functions which are used to manipulate time series, or
interact with ObsPy trace and stream objects

.. note::
    These functions have been refactored from the original SeisFlows but have
    note been tested and frankly I am not sure what the intended use is
    for the mute function is as it has not been documented or employed in
    the SeisFlows example problems.
"""
import numpy as np


def mask(slope, const, offset, nt, dt, length=400):
    """
    Constructs a tapered mask that can be applied to trace to mute early or
    late arrivals. Called by the Default preprocessing module.

    .. note::
        t_mask = slope * offset + const
        itmin = t_mask - length/2
        itmax = t_mask + length/2
        t_array = [itmin, itmax]

        offset = || s - r || is distance between source and receiver [m]
        const has units of time [s]
        slope has units of time/dist (or velocity**-1) [s/m]

    :type slope: float
    :param slope: slope applied to source receiver distance to mute arrivals
    :type const: float
    :param const: a constant time offset used to shift the mask in time
    :type offset: float
    :param offset: source-receiver distance in units of distance
    :type nt: int
    :param nt: number of samples in the waveform to be masked
    :type dt: float
    :param dt: sampling rate of the waveform to be masked
    :type length: int
    :param length: length, in time of the output mask function
    :rtype: np.array
    :return: A mask array that can be directly multipled with a waveform
    """
    # Set up the data array
    mask_arr = np.ones(nt)

    # construct taper
    win = np.sin(np.linspace(0, np.pi, 2*length))
    win = win[0:length]

    # Caculate offsets
    itmin = int(np.ceil((slope * abs(offset) + const) / dt)) - length / 2
    itmax = itmin + length

    # Generate parts of the mask array based on offsets
    if 1 < itmin < itmax < nt:
        mask_arr[0:itmin] = 0.
        mask_arr[itmin:itmax] = win * mask_arr[itmin:itmax]
    elif itmin < 1 <= itmax:
        mask_arr[0:itmax] = win[length-itmax:length] * mask_arr[0:itmax]
    elif itmin < nt < itmax:
        mask_arr[0:itmin] = 0.
        mask_arr[itmin:nt] = win[0:nt-itmin] * mask_arr[itmin:nt]
    elif itmin > nt:
        mask_arr[:] = 0.

    return mask_arr


def mute_arrivals(st, slope, const, choice):
    """
    Apply a tapered mask to a record section to mute early or late arrivals

    :type st: obspy.stream
    :param st: Stream object containing waveforms to mute
    :type slope: float
    :param slope: slope applied to source receiver distance to mute arrivals
    :type const: float
    :param const: a constant time offset used to shift the mask in time
    :type choice: str
    :param choice: "early" to mute early arrivals, "late" to mute late arrivals
    :rtype: obspy.stream
    :return: muted stream object
    """
    assert choice.upper() in ["EARLY", "LATE"]
    st_out = st.copy()

    # Get the source and receiver coordinates and time info
    nt = st_out[0].stats.npts
    dt = st_out[0].stats.delta
    s_coords = get_receiver_coords(st)
    r_coords = get_receiver_coords(st)

    for i, tr in enumerate(st_out):
        sx, sy, sz = s_coords[:][i]
        rx, ry, rz = r_coords[:][i]
        # Determine the distance between source and receiver
        offset = np.sqrt((rx - sx) ** 2 + (ry - sy) ** 2)
        mask_arr = mask(slope=slope, const=const, offset=offset, nt=nt, dt=dt)
        if choice == "early":
            tr.data *= mask_arr
        elif choice == "late":
            tr.data *= 1 - mask_arr

    return st_out


def mute_offsets(st, dist, choice):
    """
    Mute traces based on a given distance (`dist`)
    short: ||s-r|| < `dist`
    long:  ||s-r|| > `dist`

    :type st: obspy.stream
    :param st: Stream object containing waveforms to mute
    :type dist: float
    :param dist: cutoff distance
    :type choice: str
    :param choice: "short" to mute short src-rcv distances,
        "long" to mute long src-rcv distances
    :rtype: obspy.stream
    :return: muted stream object
    """
    assert choice.upper() in ["LONG", "SHORT"]
    st_out = st.copy()

    # Get the source and receiver coordinates
    s_coords = get_receiver_coords(st)
    r_coords = get_receiver_coords(st)

    for i, tr in enumerate(st_out):
        sx, sy, sz = s_coords[:][i]
        rx, ry, rz = r_coords[:][i]
        # Determine the distance between source and receiver
        offset = np.sqrt((rx - sx) ** 2 + (ry - sy) ** 2)

        if choice == "long" and (offset < dist):
            tr.data *= 0
        elif choice == "short" and (offset > dist):
            tr.data *= 0

    return st_out


def get_receiver_coords(st):
    """
    Retrieve the coordinates from a Stream object.
    Only works for SU format currently

    :type st: obspy.core.stream.Stream
    :param st: a stream to query for coordinates
    :rtype r_coords: list
    :return r_coords: list of receiver coordinates, matching the order in `st`
        ([rx], [ry], [rz])
    """
    # Seismic-Unix format
    if hasattr(st[0].stats, "su"):
        rx, ry, rz = [], [], []

        for tr in st:
            rx += [tr.stats.su.trace_header.group_coordinate_x]
            ry += [tr.stats.su.trace_header.group_coordinate_y]
            rz += [0.]
        return rx, ry, rz
    else:
        raise NotImplementedError


def get_source_coords(st):
    """
    Get the coordinates of the source object.
    Only works for SU format currently


    :type st: obspy.core.stream.Stream
    :param st: a stream to query for coordinates
    :rtype s_coords: tuple of lists
    :return s_coords: list of source coordinates, matching the order in `st`
        ([sx], [sy], [sz])
    """
    if hasattr(st[0].stats, "su"):
        sx, sy, sz = [], [], []
        for tr in st:
            sx += [tr.stats.su.trace_header.source_coordinate_x]
            sy += [tr.stats.su.trace_header.source_coordinate_y]
            sz += [0.]
        return sx, sy, sz
    else:
        raise NotImplementedError


# From the original SeisFlows code, not used but left just incase
# def correlate(u, v):
#     w = np.convolve(u, np.flipud(v))
#     return
#
#
# def tukeywin(nt, imin, imax, alpha=0.05):
#     t = np.linspace(0,1,imax-imin)
#     w = np.zeros(imax-imin)
#     p = alpha/2.
#     lo = np.floor(p*(imax-imin-1))+1
#     hi = imax-imin-lo
#     w[:lo] = (1+np.cos(np.pi/p*(t[:lo]-p)))/2
#     w[lo:hi] = np.ones((hi-lo))
#     w[hi:] = (1+np.cos(np.pi/p*(t[hi:]-p)))/2
#     win = np.zeros(nt)
#     win[imin:imax] = w
#     return win
#
# def sconvolve(s, h, w, inplace=True):
#     nt = h.nt
#     nr = h.nr
#
#     if inplace:
#         for ir in range(nr):
#             s[:,ir] = np.convolve(s[:,ir], w, 'same')
#         return s
#     else:
#         s2 = np.zeros((nt,nr))
#         for ir in range(nr):
#             s2[:,ir] = np.convolve(s[:,ir], w, 'same')
#         return s2
#
# def apply_filter_backwards(self, traces):
#     """
#     !!! This was located in seisflows.preprocess.default, but wasn't called
#     !!! anywhere. Have relinquished it to this graveyard until further notice.
#     Run the apply_filter() function but backwards
#
#     :param traces:
#     :return:
#     """
#     for tr in traces:
#         tr.data = np.flip(tr.data)
#
#     traces = self.apply_filter()
#
#     for tr in traces:
#         tr.data = np.flip(tr.data)
#
#     return traces