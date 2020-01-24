"""
SeisFlows uses obspy stream objects for holding and processing seismic data.
In some cases, obspy.read doesn't  provide the desired behavior,
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


def ascii(st, path):
    """
    Writes seismic traces as ascii files

    :type st: obspy.core.stream.Stream
    :param st: stream to write
    :type path: str
    :param path: path to datasets
    """
    for ir, tr in enumerate(st):
        nt = tr.stats.npts
        t1 = float(tr.stats.starttime)
        t2 = t1 + tr.stats.npts * tr.stats.sampling_rate

        print(nt, t1, t2)

        t = np.linspace(t1, t2, nt)
        w = tr.data

        print(os.path.join(path, tr.stats.filename))
        print(tr.times.shape, tr.data.shape)

        np.savetxt(os.path.join(path, tr.stats.filename),
                   np.column_stack((t, w))
                   )

