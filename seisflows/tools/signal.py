"""
Signal processing or manipulation functions which are used to manipulate time
series, or interact with ObsPy Trace and Stream objects. Primarily used by
the Preprocessing module
"""
import numpy as np
from seisflows import logger


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
    logger.warning("this function is currently untested, use at your own risk")

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
    Apply a tapered mask to a set of waveforms in a Stream to mute early or
    late arrivals

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
    logger.warning("this function is currently untested, use at your own risk")

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
    logger.warning("this function is currently untested, use at your own risk")

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
