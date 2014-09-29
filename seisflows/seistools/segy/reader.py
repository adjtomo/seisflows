
import numpy as np

from seisflows.tools.arraytools import uniquerows
from seisflows.tools.iotools import Reader, structure, mychar, mysize

from headers import SEGY_TAPE_LABEL, SEGY_BINARY_HEADER, SEGY_TRACE_HEADER_LONG, SEGY_TRACE_HEADER_SHORT


NMAX = 100000
FIXEDLENGTH = 1
SAVEHEADERS = 1

COORDSCALAR = 1.
DEPTHSCALAR = 1.



class SeismicReader(Reader):

  def ReadSeismicData(self):

    nsamples = int(self.read('int16',1,self.offset+114)[0])
    nbytes = int(nsamples*self.dsize+240)
    ntraces = int((self.size-self.offset)/nbytes);

    # prepare offset pointers
    if FIXEDLENGTH:
      tracelen = [nsamples]*ntraces
      traceptr = [nbytes*i + self.offset for i in range(ntraces)]

    else:
      ntraces = 1
      tracelen = []
      traceptr = []
      traceptr.append(self.offset)

      while 1:
        ntraces += 1
        nsamples = int(self.read('int16',1,traceptr[-1]+114)[0])
        nbytes = nsamples*self.dsize+240
        tracelen.append(nsamples)
        traceptr.append(traceptr[-1] + nbytes)

        if ntraces > NMAX:
          raise Exception
        elif (traceptr[-1] >= self.size):
          raise Exception
        traceptr = traceptr[:-1]
        tracelen = tracelen[:-1]

    # preallocate trace headers
    if SAVEHEADERS == 1:
      h = [self.scan(SEGY_TRACE_HEADER_LONG,traceptr[0])]
      h = h*ntraces
    elif SAVEHEADERS == 2:
      h = [self.scan(SEGY_TRACE_HEADER_SHORT,traceptr[0],contiguous=0)]
      h = h*ntraces
    else:
      h = []

    # preallocate data array
    if FIXEDLENGTH:
      d = np.zeros((nsamples,ntraces))
    else:
      d = np.zeros((tracelen.max(),len(traceptr)))

    # read trace headers and data
    for k in range(ntraces):
      if SAVEHEADERS == 1: 
        h[k] = self.scan(SEGY_TRACE_HEADER_LONG,traceptr[k])
      elif SAVEHEADERS == 2:
        h[k] = self.scan(SEGY_TRACE_HEADER_SHORT,traceptr[k],contiguous=0)

      d[:,k] = self.read(self.dtype,nsamples,traceptr[k]+240)

    self.hdrs = h
    self.data = d


  def CustomStructure(self):

    h = structure()

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

    sxyz = np.column_stack([sx,sy,sz])
    rxyz = np.column_stack([rx,ry,rz])
    ns = len(uniquerows(sxyz))
    nr = len(uniquerows(rxyz))

    # apply scaling factors
    if COORDSCALAR and DEPTHSCALAR:
      c1 = COORDSCALAR
      c2 = DEPTHSCALAR
      c3 = 1.e-6
    else:
      c1 = self.getscalar('CoordinateScalar')
      c2 = self.getscalar('ElevationOrDepthScalar')
      c3 = 1.e-6

    h.sx = c1*sx
    h.sy = c1*sy
    h.sz = c2*sz
    h.rx = c1*rx
    h.ry = c1*ry
    h.rz = c2*rz
    h.ts = c3*ts
    h.dt = c3*dt

    h.nt = nt
    h.ns = ns
    h.nr = nr

    return h


  def getarray(self,label):
    # collect array
    list = [ hdr[label] for hdr in self.hdrs ]
    return np.array(list)


  def getscalar(self,label):
    # collect scalar
    array = self.getarray(label)
    assert np.all(array == array[0])
    return array[0]




### SEGY class

class SegyReader(SeismicReader):

  def __init__(self,fname,endian=None):
    SeismicReader.__init__(self,fname,endian)

    self.dtype  = 'float'
    self.dsize  = mysize(self.dtype)
    self.offset = 0

    # check byte order
    if endian:
      self.endian = endian
    else:
      self.endian = checkByteOrder()


  def ReadSegyHeaders(self):

    # read in tape label header if present
    code = self.read('char',2,4)
    if code == 'SY':
      tapelabel = file.scan(SEGY_TAPE_LABEL,self.offset)
      self.offset += 128
    else:
      tapelabel = 'none'

    # read textual file header
    self.segyTxtHeader = self.read('char',3200,self.offset)
    self.offset += 3200

    # read binary file header
    self.segyBinHeader = self.scan(SEGY_BINARY_HEADER,self.offset)
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
      assert bool(self.segyBinHeader.FixedLengthTraceFlag)==bool(FIXEDLENGTH)



### Seismic Unix class

class SuReader(SeismicReader):

  def __init__(self,fname,endian=None):
    SeismicReader.__init__(self,fname,endian)

    self.dtype  = 'float'
    self.dsize  = mysize(self.dtype)
    self.offset = 0

    # check byte order
    if endian:
      self.endian = endian
    else:
      self.endian = checkByteOrder()



### convenience functions

def readsegy(filename):

  obj = SegyReader(filename,endian='>')

  obj.ReadSegyHeaders()
  obj.ReadSeismicData()

  d = obj.data
  h = obj.CustomStructure()

  return d, h


def readsu(filename):

  obj = SuReader(filename,endian='<')

  obj.ReadSeismicData()

  d = obj.data
  h = obj.CustomStructure()

  return d, h


