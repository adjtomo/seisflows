
from seisflows.seistools.legacy.writers import *


def su(d, path, filename):
    max_delta = 0.065535
    dummy_delta = max_delta

    if d[0].stats.delta > max_delta:
        for t in d:
            t.stats.delta = dummy_delta

    # write data to file
    d.write(path+'/'+filename, format='SU')


