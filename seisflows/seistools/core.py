
from seisflows.tools.codetools import Struct


class SeisStruct(Struct):
  """ Holds information about data
  """
  def __init__(self,nr,nt,dt,ts=0.,
               sx=[],sy=[],sz=[],
               rx=[],ry=[],rz=[],
               nrec=[],nsrc=[]):
    super(SeisStruct,self).__init__()

    self.nr = nr
    self.nt = nt
    self.dt = dt
    self.ts = ts

    self.sx = sx
    self.sy = sy
    self.sz = sz
    self.rx = rx
    self.ry = ry
    self.rz = rz

    self.nrec = nrec
    self.nsrc = nsrc


def loadascii(dir):
    wildcard = os.path.join(dir,'*.ascii')
    model = {}
    for file in glob.glob(wildcard):
      key = os.path.splitext(os.path.basename(file))[0]
      model[key] = np.loadtxt(file)
    return model


def loadsep():
    raise NotImplementedError


def loadnc():
    raise NotImplementedError


