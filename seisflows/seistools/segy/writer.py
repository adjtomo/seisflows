
import numpy as np

from seisflows.tools.iotools import Writer, mychar, mysize

from headers import SEGY_TRACE_HEADER_SHORT


class SeismicWriter(Writer):

  def __init__(self,fname):
    Writer.__init__(self,fname)

    self.dtype = 'float'
    self.dsize = mysize(self.dtype)
    self.endian = '<'
    self.offset = 0


  def prepareTraceData(self,h):

    ntraces = h.nr
    self.vals = [[1.] for k in range(ntraces)]
    c1 = 1.e0
    c2 = 1.e0
    c3 = 1.e6

    dt_max = (2.**15-1.)*10.**-6
    if h.dt >= dt_max:
      h.dt = 0

    # prepare trace headers
    for k in range(ntraces):
      self.vals[k].append(int(h.sz[k]*c1))
      self.vals[k].append(int(h.rz[k]*c1))
      self.vals[k].append(int(c1))
      self.vals[k].append(int(c2))
      self.vals[k].append(int(h.sx[k]*c2))
      self.vals[k].append(int(h.sy[k]*c2))
      self.vals[k].append(int(h.rx[k]*c2))
      self.vals[k].append(int(h.ry[k]*c2))
      self.vals[k].append(int(h.ts*c3))
      self.vals[k].append(int(h.nt))
      self.vals[k].append(int(h.dt*c3))


  def writeTraceData(self,d):

    nsamples = d.shape[0]
    nbytes = nsamples*self.dsize+240
    ntraces = d.shape[1]

    for k in range(ntraces):

      # write trace header
      self.printf(SEGY_TRACE_HEADER_SHORT,self.vals[k],k*nbytes,contiguous=0)

      # write trace data
      self.write(self.dtype,d[:,k],nsamples,k*nbytes+240)


class SuWriter(SeismicWriter):

  def __init__(self,fname):
    SeismicWriter.__init__(self,fname)


def writesegy():
  pass


def writesu(filename,d,h):

  obj = SuWriter(filename)

  obj.prepareTraceData(h)
  obj.writeTraceData(d)

