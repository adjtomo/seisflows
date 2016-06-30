import numpy as np

from seisflows.tools.io import BinaryWriter, mychar, mysize

from headers import SEGY_TRACE_HEADER

FIELDS = [
    'TraceSequenceLine',
    'SourceWaterDepth',
    'GroupWaterDepth',
    'ElevationOrDepthScalar',
    'CoordinateScalar',
    'SourceX',
    'SourceY',
    'GroupX',
    'GroupY',
    'RecordingDelay_ms',
    'NumberSamples',
    'SampleInterval_ms']

# cull header fields
_tmp = []
for field in SEGY_TRACE_HEADER:
    if field[-1] in FIELDS:
        _tmp.append(field)
SEGY_TRACE_HEADER = _tmp


class SeismicWriter(BinaryWriter):
    def __init__(self, fname):
        super(SeismicWriter, self).__init__(fname)

        self.dtype = 'float'
        self.dsize = mysize(self.dtype)
        self.endian = '<'
        self.offset = 0

    def prepareTraceData(self, h):

        nr = int(h.nr)
        self.vals = [[1] for k in range(nr)]
        self.ntraces = nr

        c1 = 1
        c2 = 1
        c3 = 1000000

        nt = int(h.nt)
        dt = int(h.dt*c3)
        ts = int(h.ts*c3)

        dt_max = (2.**15 - 1.)*10.**-6
        if h.dt >= dt_max:
            dt = 0

        sx = self.getarray(h, 'sx', c1)
        sy = self.getarray(h, 'sy', c1)
        sz = self.getarray(h, 'sz', c2)

        rx = self.getarray(h, 'rx', c1)
        ry = self.getarray(h, 'ry', c1)
        rz = self.getarray(h, 'rz', c2)

        # prepare trace headers
        for k in range(nr):
            self.vals[k].append(sz[k])
            self.vals[k].append(rz[k])
            self.vals[k].append(c1)
            self.vals[k].append(c2)
            self.vals[k].append(sx[k])
            self.vals[k].append(sy[k])
            self.vals[k].append(rx[k])
            self.vals[k].append(ry[k])
            self.vals[k].append(ts)
            self.vals[k].append(nt)
            self.vals[k].append(dt)

    def getarray(self, h, key, constant):
        if h[key] == []:
            return [0]*self.ntraces

        array = [int(f*constant) for f in h[key]]
        return array

    def writeTraceData(self, d):

        nsamples = d.shape[0]
        nbytes = nsamples*self.dsize + 240
        nr = d.shape[1]

        for k in range(nr):
            # write trace header
            self.printf(SEGY_TRACE_HEADER, self.vals[k],
                        k*nbytes, contiguous=False)

            # write trace data
            self.write(self.dtype, d[:, k], nsamples, k*nbytes + 240)


class SuWriter(SeismicWriter):
    def __init__(self, fname):
        SeismicWriter.__init__(self, fname)


def writesegy():
    raise NotImplementedError


def writesu(filename, d, h):
    obj = SuWriter(filename)
    obj.prepareTraceData(h)
    obj.writeTraceData(d)
