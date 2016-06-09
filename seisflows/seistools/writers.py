
from seisflows.seistools.legacy.writers import *


def su(d, path, filename):
    max_delta = 0.065535
    dummy_delta = max_delta

    if d[0].stats.delta > max_delta:
        for t in d:
            t.stats.delta = dummy_delta

    # write data to file
    d.write(path+'/'+filename, format='SU')


def ascii(stream, path, filenames):
    import numpy as np

    for tr in stream:
        t1 = tr.starttime
        t2 = tr.starttime + tr.npts*tr.sampling_rate
        times = np.arange(t1, t2, tr.npts)
        np.savetxt(path +'/'+ tr.filename,
                   np.column_stack(times, tr.data))
