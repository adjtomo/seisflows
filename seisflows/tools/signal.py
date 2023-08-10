"""
Signal processing or manipulation functions which are used to manipulate time
series, or interact with ObsPy Trace and Stream objects. Primarily used by
the Preprocessing module
"""
import numpy as np
from seisflows import logger


def filter(st, choice, min_freq=None, max_freq=None, zerophase=True, **kwargs):
    """
    Apply a filter to waveform data using ObsPy, throw on a standard
    demean, detrened and taper prior to filtering. Options for different
    filtering types. Uses default filter options from ObsPy.

    Zerophase enforced to be True to avoid phase shifting data.

    :type st: obspy.core.stream.Stream
    :param st: stream to be filtered
    :rtype: obspy.core.stream.Stream
    :return: filtered traces
    """
    if choice.upper() == "BANDPASS":
        st.filter("bandpass", zerophase=zerophase,
                  freqmin=min_freq, freqmax=max_freq)
    elif choice.upper() == "LOWPASS":
        st.filter("lowpass", zerophase=zerophase, freq=max_freq)
    elif choice.upper() == "HIGHPASS":
        st.filter("highpass", zerophase=zerophase, freq=min_freq)
    else:
        raise NotImplementedError(f"filter choice {choice} is not available")

    return st


def resample(st_a, st_b):
    """
    Resample all traces in `st_a` to the sampling rate of `st_b`. Resamples
    one to one, that is each trace in `obs` is resampled to the
    corresponding indexed trace in `syn`

    :type st_a: obspy.core.stream.Stream
    :param st_a: stream to be resampled using sampling rates from `st_b`
    :type st_b: obspy.core.stream.Stream
    :param st_b:  stream whose sampling rates will be used to resample
        `st_a`. Usually this is the synthetic data
    :rtype: (obspy.core.stream.Stream, obspy.core.stream.Stream)
    :return: `st_a` (resampled), `st_b` (original)
    """
    for tr_a, tr_b in zip(st_a, st_b):
        sr_a = tr_a.stats.sampling_rate
        sr_b = tr_b.stats.sampling_rate
        if sr_a != sr_b:
            logger.debug(f"resampling '{tr_a.get_id()}' {sr_a}->{sr_b} Hz")
            tr_a.resample(sampling_rate=tr_b.stats.sampling_rate)

    return st_a, st_b


def normalize(st, choice=None, st_rel=None):
    """
    Normalize amplitudes of Stream object waveforms based on the given choice of
    normalization function.

    :type st: obspy.core.stream.Stream
    :param st: All of the data streams to be normalized
    :type choice: str
    :param choice: choice of normalization parameter, from the following:
        - None: Do not normalize. Used to bypass procedure
        TRACE-WISE NORMALIZATION
        - TNORML1: normalize per trace by the L1 norm of itself
        - TNORML2: normalize per trace by the L2 norm of itself
        - TNORM_MAX: normalize by the maximum positive amplitude in the trace
        - TNORM_ABSMAX: normalize by the absolute maximum amplitude in the trace
        - TNORM_MEAN: normalize by the mean of the absolute trace
        RELATIVE NORMALIZATION
        - RNORM_MAX: normalize `st` by the max positive amplitude of `st_rel`
        - RNORM_ABSMAX: normalize `st` by abs max amplitude of `st_rel`
    :type st_rel: obspy.core.stream.Stream
    :param st_rel: Second stream used for relative normalization. Optional
        and only required if 'RNORM' set as `choice`
    :rtype: obspy.core.stream.Stream
    :return: stream with normalized traces
    """
    if choice is None:
        return

    choice = choice.upper()
    st_out = st.copy()
    if "RNORM" in choice:
        assert(st_rel is not None), (
            f"corresponding 'relative' stream is required for RNORM type "
            f"normalizations"
        )

    # Normalize each trace by its L1 norm
    if normalize == "TNORML1":
        for tr in st_out:
            w = np.linalg.norm(tr.data, ord=1)
            if w < 0:
                logger.warning(f"CAUTION: L1 Norm for {tr.get_id()} is "
                               f"negative, this will result in "
                               f"unintentional sign flip")
            tr.data /= w
    # Normalize each trace by its L2 norm
    elif normalize == "TNORML2":
        for tr in st_out:
            w = np.linalg.norm(tr.data, ord=2)
            if w < 0:
                logger.warning(f"CAUTION: L2 Norm for {tr.get_id()} is "
                               f"negative, this will result in "
                               f"unintentional sign flip")
            tr.data /= w
    # Normalize each trace by its maximum positive amplitude
    elif normalize == "TNORM_MAX":
        for tr in st_out:
            w = np.max(tr.data)
            tr.data /= w
    # Normalize each trace by the maximum amplitude (neg or pos)
    elif normalize == "TNORM_ABSMAX":
        for tr in st_out:
            w = np.abs(tr.max())
            tr.data /= w
    # Normalize by the mean of absolute trace amplitudes
    elif normalize == "TNORM_MEAN":
        for tr in st_out:
            w = np.mean(np.abs(tr.data))
            tr.data /= w
    # Normalize by the max amplitude of the corresponding relative trace
    elif normalize == "RNORM_MAX":
        for tr, tr_rel in zip(st_out, st_rel):
            w = np.max(tr_rel)
            tr.data /= w
    elif normalize == "RNORM_MAX":
        for tr, tr_rel in zip(st_out, st_rel):
            w = np.abs(np.max(tr_rel))
            tr.data /= w
    else:
        raise NotImplementedError(f"normalization choice '{choice}' is not "
                                  f"implemented")

    # !!! These are not currently working. Open a GitHub issue if you
    # !!! would like to see event-wise normalization
    # Normalize an event by the L1 norm of all traces
    # if 'ENORML1' in norm_choices:
    #     w = 0.
    #     for tr in st_out:
    #         w += np.linalg.norm(tr.data, ord=1)
    #     for tr in st_out:
    #         tr.data /= w
    # # Normalize an event by the L2 norm of all traces
    # elif "ENORML2" in norm_choices:
    #     w = 0.
    #     for tr in st_out:
    #         w += np.linalg.norm(tr.data, ord=2)
    #     for tr in st_out:
    #         tr.data /= w

    return st_out


def mute():
    """
    Mute arrivals or offsets on a waveform based on slopes and constants
    defined by the User.
    """


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

        for tr in st_out:
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
        for tr in st_out:
            sx += [tr.stats.su.trace_header.source_coordinate_x]
            sy += [tr.stats.su.trace_header.source_coordinate_y]
            sz += [0.]
        return sx, sy, sz
    else:
        raise NotImplementedError
