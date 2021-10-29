"""
SeisFlows uses obspy stream objects for holding and processing seismic data.
In some cases, obspy.read doesn't provide the desired behavior,
so we introduce an additonal level of indirection

Used by the PREPROCESS class and specified by the WRITER parameter
"""
import os
import numpy as np


def su(st, path, filename):
    """
    Writes seismic unix files outputted by Specfem, using Obspy

    :type st: obspy.core.stream.Stream
    :param st: stream to write
    :type path: str
    :param path: path to datasets
    :type filename: str
    :param filename: file to read
    """
    for tr in st:
        # Work around obspy data type conversion
        tr.data = tr.data.astype(np.float32)

    max_delta = 0.065535
    dummy_delta = max_delta

    if st[0].stats.delta > max_delta:
        for tr in st:
            tr.stats.delta = dummy_delta

    # Write data to file
    st.write(os.path.join(path, filename), format='SU')


def ascii(st, path, filename=None):
    """
    Writes seismic traces as ascii files. Kwargs are left to keep structure of 
    inputs compatible with other input formats.

    :type st: obspy.core.stream.Stream
    :param st: stream to write
    :type path: str
    :param path: path to datasets
    """
    for tr in st:
        if filename is None:
            filename = tr.stats.filename

        fid_out = os.path.join(path, filename)
        
        # Float provides the time difference between starttime and default time
        time_offset = float(tr.stats.starttime)

        data_out = np.vstack((tr.times() + time_offset, tr.data)).T

        np.savetxt(fid_out, data_out, ["%13.7f", "%17.7f"])

