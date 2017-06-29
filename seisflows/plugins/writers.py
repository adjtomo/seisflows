
# SeisFlows uses obspy stream objects for holding and processing seismic data. In some cases, obspy.read doesn't  provide the desired behavior, so we introduce an additonal level of indirection

# used by the PREPROCESS class and specified by the WRITER parameter


import numpy as np

def su(d, path, filename):
    for t in d:
        # work around obspy data type conversion
        t.data = t.data.astype(np.float32)


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

        t = np.linspace(t1, t2, nt)
        w = tr.data

        print path +'/'+ tr.stats.filename
        print times.shape, tr.data.shape
        np.savetxt(path +'/'+ tr.stats.filename,
                   np.column_stack((t, w)))

