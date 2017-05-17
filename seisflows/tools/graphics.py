
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from obspy.core.stream import Stream


def plot_gll(x, y, z):
    """ Plots values on 2D unstructured GLL mesh
    """
    r = (max(x) - min(x))/(max(y) - min(y))
    rx = r/np.sqrt(1 + r**2)
    ry = 1/np.sqrt(1 + r**2)

    f = plt.figure(figsize=(10*rx, 10*ry))
    p = plt.tricontourf(x, y, z, 125)
    plt.axis('image')
    return f, p


def plot_vector(t, v, xlabel='', ylabel='', title=''):
    """ Plots a vector or time series.

    Parameters
    ----------
    v: ndarray, ndims = 1/2
        Vector or time series to plot
    xlabel: str
        x axis label
    ylabel: str
        y axis label
    title: str
        plot title

    Raises
    ------
    ValueError
        If dimensions of v are greater than 2
    """

    # check input dimension
    if v.ndim > 2:
        raise ValueError('v must be a vector or a time series')

    if v.ndim == 1:
        x = range(len(v))
        y = v
    else:
        x = v[:, 0]
        y = v[:, 1]

    # plot
    plt.plot(t, v)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()


def plot_section(stream, ax=None, cmap='seismic', clip=100, title='', x_interval=1.0, y_interval=1.0):
    """  Plots a seismic section from an obspy stream.

    Parameters
    ----------
    stream: Obspy stream object
        Obspy stream object created from a SU data file
    ax: Matplotlib Axes object
        Optional axis object
    cmap: str
        Matplotlib colormap option.
    clip: float
        Percentage value (0-100) for amplitude clipping
    title: str
        plot title
    x_interval: float
        Offset axis tick interval in km
    y_interval: float
        Time axis tick interval in km

    Raises
    ------
    NotImplementedError
        If stream object does not have SU format
    """

    # check format of stream
    if stream[0].stats._format != 'SU':
        raise NotImplemented('plot_section currently only supports streams for SU data files.')

    # get dimensions
    nr = len(stream)
    nt = len(stream[0].data)
    dt = stream[0].stats.delta
    d_aspect = nr / float(nt)

    # convert stream to image array
    data = _convert_to_array(stream)

    # default values
    fsize = 6
    scale_factor = 1.5

    if ax is None:
        fig, ax = plt.subplots(figsize=(fsize, scale_factor*fsize))

    im = ax.imshow(data, aspect=scale_factor*d_aspect, clim=_cscale(data, clip=clip))
    im.set_cmap(cmap)

    # labels
    ax.set_title(title)
    ax.set_xlabel('Offset [km]')
    ax.set_ylabel('Time [s]')

    #set ticks
    t = _get_time(stream)
    yticks, ytick_labels = get_regular_ticks(t, y_interval)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels)

    offsets =_get_offsets(stream)
    xticks, xtick_labels = get_regular_ticks(offsets, x_interval)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)

    return ax


def _convert_to_array(stream):
    """ Extracts trace data from an obspy stream and returns a 2D array.

    Parameters
    ----------
    stream: Obspy stream object
        Stream storing trace data

    Returns
    -------
    output: ndarray, ndim=2
        Returns an (nt*nr) array. nt and nr are the number of sample points
        and number of traces respectively. Assumes trace lengths are equal
        for all traces.

    Raises
    ------
    TypeError
        If stream is not an obspy stream
    """
    if not isinstance(stream, Stream):
        raise TypeError('Input object should be an obspy stream.')

    nt = len(stream.traces[0].data)
    nr = len(stream)
    output = np.zeros((nt, nr))

    for i, trace in enumerate(stream):
        output[:, i] = trace.data[:]

    return output


def _cscale(v, clip=100):
    """ Return limits for colormap.
    """
    perc = clip / 100.
    return -perc * abs(v).max(), perc * abs(v).max()


def _get_time(stream):
    """ Get fixed time vector for stream object.
    """
    dt = stream[0].stats.delta
    nt = len(stream[0].data)
    return np.arange(0, nt*dt, dt)


def _get_offsets(stream):
    """ Return offsets.
    """
    nr = len(stream)
    offsets = np.zeros(nr)
    scalco = stream[0].stats.su.trace_header.scalar_to_be_applied_to_all_coordinates

    # set scale to km
    if scalco == 0:
        scalco = 1e-3 # assume coords are in m
    else:
        scalco = 1.0e-3 / scalco

    for i, tr in enumerate(stream):
        offsets[i] = (tr.stats.su.trace_header.group_coordinate_x -
                      tr.stats.su.trace_header.source_coordinate_x) * scalco
    return offsets


def get_regular_ticks(v, interval):
    """ Returns regular tick intervals.
    """
    f = interp1d(v, range(len(v)))
    begin = int(v[0] / interval) * interval
    end = v[-1]
    tick_labels = np.arange(begin, end, interval)
    ticks = f(tick_labels)

    return ticks, tick_labels
