
from collections import Mapping

from seisflows.tools.code import Struct


class SeisStruct(Struct):
    """ Holds information about data
    """

    def __init__(self, nr=0, nt=0, dt=0., ts=0.,
                 sx=[], sy=[], sz=[],
                 rx=[], ry=[], rz=[],
                 nrec=[], nsrc=[]):
        super(SeisStruct, self).__init__(
            [['nr', nr], ['nt', nt], ['dt', dt], ['ts', ts],
             ['sx', sx], ['sy', sy], ['sz', sz],
             ['rx', rx], ['ry', ry], ['rz', rz],
             ['nrec', nrec], ['nsrc', nsrc]])


class ModelStruct(Mapping):
    def __init__(self, reader):
        self.reader = reader
        self.cached_key = None
        self.cached_val = None

    def __getitem__(self, key):
        if key == self.cached_key:
            return self.cached_val
        else:
            array = self.reader(key)
            self.cached_key = key
            self.cached_val = array
        return array

    def __iter__(self):
        return []

    def __len__(self):
        return len([])


def loadascii(dir):
    wildcard = os.path.join(dir, '*.ascii')
    model = {}
    for file in glob.glob(wildcard):
        key = os.path.splitext(os.path.basename(file))[0]
        model[key] = np.loadtxt(file)
    return model


def loadsep():
    raise NotImplementedError


def loadnc():
    raise NotImplementedError
