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


def su(filename):
    """
    Reads seismic unix files outputted by Specfem, using Obspy

    :type filename: str
    :param filename: full path to data file to read
    """
    st = read(os.path.join(filename), format='SU', byteorder='<')
    
    return st


def ascii(filename):
    """
    Reads SPECFEM3D-style ASCII data

    :type filename: str
    :param filename: full path to data file to read
    """
    st = Stream()
    stats = Stats()

    time, data = loadtxt(filename).T

    stats.filename = os.path.basename(filename)
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

    st.append(Trace(data=data, header=stats))

    return st

