
from seisflows.seistools.legacy.readers import *


def su(path, filename):
    import obspy
    s = obspy.read(path +'/'+ filename, 
                   format='SU',
                   byteorder='<')
    return s


def segy(path, filename):
    import obspy
    s = obspy.read(path +'/'+ filename,   
                   format='SEGY',
                   byteorder='<')
    return s


