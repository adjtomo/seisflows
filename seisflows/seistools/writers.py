
import numpy as np

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
    for ir, tr in enumerate(stream):
        nt = tr.stats.npts
        t1 = float(tr.stats.starttime)
        t2 = t1 + tr.stats.npts*tr.stats.sampling_rate
        print nt, t1, t2
        times = np.linspace(t1, t2, nt)
        print path +'/'+ tr.stats.filename
        print times.shape, tr.data.shape
        np.savetxt(path +'/'+ tr.stats.filename,
                   np.column_stack((times, tr.data)))

