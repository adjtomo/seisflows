
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
    def __init__(self):
        raise NotImplementedError


def load1(dirname, parameters, mapping, nproc, logfile=None):
    """ reads SPECFEM model

      Models are stored as a Fortran binary fomrat and and separated into 
      mulitple files according to material parameter and processor rank.
    """
    parts = {}
    for key in sorted(parameters):
        parts[key] = []
        for iproc in range(nproc):
            filename = 'proc%06d_%s.bin' % (iproc, mapping(key))
            part = loadbin(join(dirname, filename))
            parts[key].append(part)
    return parts


def load2(dirname, parameters, mapping, nproc, logfile):
    """ reads SPECFEM model

      Provides the same functionality as load1 but with debugging output.
    """
    # read database files
    parts = {}
    minmax = {}
    for key in sorted(parameters):
        parts[key] = []
        minmax[mapping(key)] = [+np.Inf,-np.Inf]
        for iproc in range(nproc):
            filename = 'proc%06d_%s.bin' % (iproc, mapping(key))
            part = loadbin(join(dirname, filename))
            parts[key].append(part)
            # keep track of min, max
            if part.min() < minmax[mapping(key)][0]: minmax[mapping(key)][0] = part.min()
            if part.max() > minmax[mapping(key)][1]: minmax[mapping(key)][1] = part.max()

    # print min, max
    if logfile:
        with open(logfile,'a') as f:
            f.write(abspath(dirname)+'\n')
            for key,val in minmax.items():
                f.write('%-15s %10.3e %10.3e\n' % (key, val[0], val[1]))
            f.write('\n')
    return parts


def load3(dirname, parameters, mapping, nproc, logfile):
    """ reads SPECFEM model

      Provides the same functionality as load1 but with debugging output and
      improved memory usage.
    """
    minmax = {}

    def helper_func(key):
        parts = []
        minmax[mapping(key)] = [+np.Inf,-np.Inf]
        for iproc in range(nproc):
            filename = 'proc%06d_%s.bin' % (iproc, mapping(key))
            part = loadbin(join(dirname, filename))
            parts.append(part)
            # keep track of min, max
            if part.min() < minmax[mapping(key)][0]: minmax[mapping(key)][0] = part.min()
            if part.max() > minmax[mapping(key)][1]: minmax[mapping(key)][1] = part.max()
        return parts

    class helper_class(Mapping):
        def __init__(self):
            self.cached_key = None
            self.cached_val = None

        def __getitem__(self, key):
            if key == self.cached_key:
                return self.cached_val
            else:
                parts = helper_func(key)
                self.cached_key = key
                self.cached_val = parts
            return parts

        def __iter__(self):
            return []

        def __len__(self):
            return len([])

        def __del__(self):
            if logfile:
                with open(logfile,'a') as f:
                    f.write(abspath(dirname)+'\n')
                    for key in sorted(minmax):
                        val = minmax[key]
                        f.write('%-15s %10.3e %10.3e\n' % (key, val[0], val[1]))
                    f.write('\n')

    return helper_class()


load = load3

