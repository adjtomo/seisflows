
# SeisFlows uses obspy stream objects for holding and processing seismic data. In some cases, obspy.read doesn't  provide the desired behavior, so we introduce an additonal level of indirection

# used by the PREPROCESS class and specified by the READER parameter


def su(path, filename):
    """ Reads Seismic Unix files

      Hardwired ''
    """
    import obspy
    stream = obspy.read(path +'/'+ filename, 
                   format='SU',
                   byteorder='<')
    return stream


def ascii(path, filenames):
    """ Reads SPECFEM3D-style ascii data
    """
    from numpy import loadtxt
    from obspy.core import Stream, Stats, Trace

    stream = Stream()
    for filename in filenames:
        stats = Stats()
        data = loadtxt(path +'/'+ filename)

        stats.filename = filename
        stats.starttime = data[0,0]
        stats.sampling_rate = data[0,1] - data[0,0]
        stats.npts = len(data[:,0])

        try:
            parts = filename.split('.')
            stats.network = parts[0]
            stats.station = parts[1]
            stats.channel = temp[2]
        except:
            pass

        stream.append(Trace(data=data[:,1], header=stats))

    return stream

