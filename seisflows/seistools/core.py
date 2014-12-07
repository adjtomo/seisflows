
from seisflows.tools.codetools import Struct


class SeisStruct(Struct):
    """ Holds information about data
    """
    def __init__(self,nr=0,nt=0,dt=0.,ts=0.,
                 sx=[],sy=[],sz=[],
                 rx=[],ry=[],rz=[],
                 nrec=[],nsrc=[]):

        super(SeisStruct,self).__init__([['nr',nr],['nt',nt],['dt',dt],['ts',ts],
                 ['sx',sx],['sy',sy],['sz',sz],
                 ['rx',rx],['ry',ry],['rz',rz],
                 ['nrec',nrec],['nsrc',nsrc]])

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
