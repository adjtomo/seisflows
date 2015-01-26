import numpy as np

from seisflows.tools.array import uniquerows
from seisflows.tools.code import Struct
from seisflows.tools.io import BinaryReader, mychar, mysize

from seisflows.seistools.shared import SeisStruct
from seisflows.seistools.segy.headers import \
    SEGY_TAPE_LABEL, SEGY_BINARY_HEADER, SEGY_TRACE_HEADER


NMAX = 100000
FIXEDLENGTH = True
SAVEHEADERS = True
COORDSCALAR = 1.
DEPTHSCALAR = 1.

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


class SeismicReader(BinaryReader):
    """ Base class used by both SegyReader and SuReader
    """

    def ReadSeismicData(self):
        nsamples = int(self.read('int16', 1, self.offset + 114)[0])
        nbytes = int(nsamples*self.dsize + 240)
        ntraces = int((self.size - self.offset)/nbytes)

        # prepare offset pointers
        if FIXEDLENGTH:
            tracelen = [nsamples]*ntraces
            traceptr = [nbytes*i + self.offset for i in range(ntraces)]

        else:
            ntraces = 1
            tracelen = []
            traceptr = [self.offset]

            while 1:
                ntraces += 1
                nsamples = int(self.read('int16', 1, traceptr[-1] + 114)[0])
                nbytes = nsamples*self.dsize + 240
                tracelen.append(nsamples)
                traceptr.append(traceptr[-1] + nbytes)

                if ntraces > NMAX:
                    raise Exception
                elif traceptr[-1] >= self.size:
                    raise Exception
                traceptr = traceptr[:-1]
                tracelen = tracelen[:-1]

        # preallocate trace headers
        if SAVEHEADERS:
            h = [self.scan(SEGY_TRACE_HEADER, traceptr[0], contiguous=False)]
            h = h*ntraces
        else:
            h = []

        # preallocate data array
        if FIXEDLENGTH:
            d = np.zeros((nsamples, ntraces))
        else:
            d = np.zeros((tracelen.max(), len(traceptr)))

        # read trace headers and data
        for k in range(ntraces):
            if SAVEHEADERS:
                h[k] = self.scan(SEGY_TRACE_HEADER, traceptr[k],
                                 contiguous=False)
            d[:, k] = self.read(self.dtype, nsamples, traceptr[k] + 240)

        # store results
        self.ntraces = ntraces
        self.hdrs = h
        self.data = d

    def getstruct(self):
        nr = self.ntraces

        # collect scalars
        nt = self.getscalar('NumberSamples')
        ts = self.getscalar('RecordingDelay_ms')
        dt = self.getscalar('SampleInterval_ms')

        # collect arrays
        sx = self.getarray('SourceX')
        sy = self.getarray('SourceY')
        sz = self.getarray('SourceWaterDepth')

        rx = self.getarray('GroupX')
        ry = self.getarray('GroupY')
        rz = self.getarray('GroupWaterDepth')

        # apply scaling factors
        if COORDSCALAR and DEPTHSCALAR:
            c1 = COORDSCALAR
            c2 = DEPTHSCALAR
            c3 = 1.e-6
        else:
            c1 = self.getscalar('CoordinateScalar')
            c2 = self.getscalar('ElevationOrDepthScalar')
            c3 = 1.e-6

        sxyz = np.column_stack([sx, sy, sz])
        rxyz = np.column_stack([rx, ry, rz])
        nsrc = len(uniquerows(sxyz))
        nrec = len(uniquerows(rxyz))

        return SeisStruct(nr, nt, dt, ts,
                          c1*sx, c1*sy, c2*sz,
                          c1*rx, c1*ry, c2*rz,
                          nsrc, nrec)

    def getarray(self, key):
        # collect array
        list = [hdr[key] for hdr in self.hdrs]
        return np.array(list)

    def getscalar(self, key):
        # collect scalar
        array = self.getarray(key)
        return array[0]


class SegyReader(SeismicReader):
    """ SEGY reader
    """

    def __init__(self, fname, endian=None):
        SeismicReader.__init__(self, fname, endian)

        self.dtype = 'float'
        self.dsize = mysize(self.dtype)
        self.offset = 0

        # check byte order
        if endian:
            self.endian = endian
        else:
            raise ValueError("SU Reader should specify the endianness")

    def ReadSegyHeaders(self):
        # read in tape label header if present
        code = self.read('char', 2, 4)
        if code == 'SY':
            tapelabel = file.scan(SEGY_TAPE_LABEL, self.offset)
            self.offset += 128
        else:
            tapelabel = 'none'

        # read textual file header
        self.segyTxtHeader = self.read('char', 3200, self.offset)
        self.offset += 3200

        # read binary file header
        self.segyBinHeader = self.scan(SEGY_BINARY_HEADER, self.offset)
        self.offset += 400

        # read in extended textual headers if present

        self.CheckSegyHeaders()

    def CheckSegyHeaders(self):
        # check revision number
        self.segyvers = '1.0'

        # check format code
        self.segycode = 5

        # check trace length
        if FIXEDLENGTH:
            assert bool(self.segyBinHeader.FixedLengthTraceFlag) == bool(
                FIXEDLENGTH)


class SuReader(SeismicReader):
    """ Seismic Unix file reader
    """

    def __init__(self, fname, endian=None):
        SeismicReader.__init__(self, fname, endian)

        self.dtype = 'float'
        self.dsize = mysize(self.dtype)
        self.offset = 0

        # check byte order
        if endian:
            self.endian = endian
        else:
            raise ValueError("SU Reader should specify the endianness")


def readsegy(filename):
    """ SEGY convenience function
    """
    obj = SegyReader(filename, endian='>')
    obj.ReadSegyHeaders()
    obj.ReadSeismicData()

    d = obj.data
    h = obj.getstruct()
    return d, h


def readsu(filename):
    """ SU convenience function
    """
    obj = SuReader(filename, endian='<')
    obj.ReadSeismicData()

    d = obj.data
    h = obj.getstruct()
    return d, h
