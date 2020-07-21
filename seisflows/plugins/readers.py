"""
SeisFlows uses obspy stream objects for holding and processing seismic data.
In some cases, obspy.read doesn't provide the desired behavior, so we
introduce an additonal level of indirection

Used by the PREPROCESS class and specified by the READER parameter
"""
import os


def su(path, filename):
    """
    Reads seismic unix files outputted by Specfem, using Obspy

    :type path: str
    :param path: path to datasets
    :type filename: str
    :param filename: file to read
    """
    from obspy import read

    st = read(os.path.join(path, filename), format='SU', byteorder='<')
    
    return st


def ascii(path, filenames):
    """
    Reads SPECFEM3D-style ASCII data

    :type path: str
    :param path: path to datasets
    :type filenames: list
    :param filenames: files to read
    """
    from numpy import loadtxt
    from obspy.core import Stream, Stats, Trace

    stream = Stream()
    for filename in filenames:
        stats = Stats()
        data = loadtxt(os.path.join(path, filename))

        stats.filename = filename
        stats.starttime = data[0, 0]
        stats.sampling_rate = data[0, 1] - data[0, 0]
        stats.npts = len(data[:, 0])

        try:
            parts = filename.split('.')
            stats.network = parts[0]
            stats.station = parts[1]
            stats.channel = parts[2]
        except:
            pass

        stream.append(Trace(data=data[:, 1], header=stats))

    return stream

