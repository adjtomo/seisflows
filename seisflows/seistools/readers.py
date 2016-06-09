
from seisflows.seistools.legacy.readers import *


def su(path, filename):
    import obspy
    stream = obspy.read(path +'/'+ filename, 
                   format='SU',
                   byteorder='<')
    return stream


def segy(path, filename):
    import obspy
    stream = obspy.read(path +'/'+ filename,   
                   format='SEGY',
                   byteorder='<')
    return stream


def ascii(path, filenames):
    from numpy import loadtxt
    from obspy.core import Stream, Stats, Trace

    stream = Stream()
    for filename in filenames:
        stats = Stats()
        data = loadtxt(path +'/'+ filename)

        stats.filename = filename
        stats.starttime = data[0,0]
        stats.sampling_rate = data[0,0] - data[0,1]
        stats.npts = len(data[0,:])

        #stats.network = temp[0]
        #stats.station = temp[1]
        #stats.location = temp[2]
        #stats.channel = temp[3]

        stream.append(Trace(data=data[:,1], header=stats))

    return stream

