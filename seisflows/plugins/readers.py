"""
SeisFlows uses obspy stream objects for holding and processing seismic data.
In some cases, obspy.read doesn't provide the desired behavior, so we
introduce an additonal level of indirection

Used by the PREPROCESS class and specified by the READER parameter
"""
import os
from numpy import loadtxt
from obspy import read 
from obspy.core import Stream, Stats, Trace


def su(path, filename):
    """
    Reads seismic unix files outputted by Specfem, using Obspy

    :type path: str
    :param path: path to datasets
    :type filename: str
    :param filename: file to read
    """
    st = read(os.path.join(path, filename), format='SU', byteorder='<')
    
    return st


def ascii(path, filename):
    """
    Reads SPECFEM3D-style ASCII data

    :type path: str
    :param path: path to datasets
    :type filenames: list
    :param filenames: files to read
    """
    stream = Stream()
    stats = Stats()

    time, data = loadtxt(os.path.join(path, filename)).T

    stats.filename = filename
    stats.starttime = time[0]
    stats.delta = time[1] - time[0]
    stats.npts = len(data)

    try:
        parts = filename.split(".")
        stats.network = parts[0]
        stats.station = parts[1]
        stats.channel = parts[2]
    except:
        pass

    stream.append(Trace(data=data, header=stats))

    return stream

