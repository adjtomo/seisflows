
import numpy as np

from collections import Mapping

from seisflows.tools.code import Struct, abspath, join
from seisflows.tools.io import loadbin


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


def load(dirname, parameters, mapping, nproc, output):
    """ reads SPECFEM model
    """
    # read database files
    parts = {}
    minmax = {}
    for key in parameters:
        parts[key] = []
        minmax[key] = [+np.Inf,-np.Inf]
        for iproc in range(nproc):
            filename = 'proc%06d_%s.bin' % (iproc, mapping(key))
            part = loadbin(join(dirname, filename))
            parts[key].append(part)
            # keep track of min, max
            if part.min() < minmax[key][0]: minmax[key][0] = part.min()
            if part.max() > minmax[key][1]: minmax[key][1] = part.max()

    # print min, max
    if output:
        with open(output,'a') as f:
            f.write(abspath(dirname)+'\n')
            for key,val in minmax.items():
                f.write('%-10s %10.3e %10.3e\n' % (key, val[0], val[1]))
            f.write('\n')
    return parts

