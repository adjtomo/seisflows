
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
        self.keys = set()

    def __getitem__(self, key, verbose=False):
        self.keys.add(key)
        val = self.reader(key)

        if verbose:
            print '%s: %e %e' % (key, val.min(), val.max())

        return val

    def __iter__(self):
        return self.keys

    def __len__(self):
        return len(self.keys)


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
