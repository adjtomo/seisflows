
from string import find

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
    def __init__(self):
        raise NotImplementedError


def getpar(key, file='DATA/Par_file', sep='='):
    """ Reads parameter from SPECFEM parfile
    """
    with open(file, 'r') as f:

        # read line by line
        for line in f:
            if find(line, key) == 0:
                # read key
                key, val = _split(line, sep)
                if not key:
                    continue
                # read val
                val, _ = _split(val, '#')
                return val.strip()

        raise Exception


def setpar(key, val, file='DATA/Par_file', path='.', sep='='):
    """ Writes parameter to SPECFEM parfile
    """

    val = str(val)

    # read line by line
    with open(path + '/' + file, 'r') as f:
        lines = []
        for line in f:
            if find(line, key) == 0:
                # read key
                key, _ = _split(line, sep)
                # read comment
                _, comment = _split(line, '#')
                n = len(line) - len(key) - len(val) - len(comment) - 2
                # replace line
                if comment:
                    line = _merge(key, sep, val, ' '*n, '#', comment)
                else:
                    line = _merge(key, sep, str(val), '\n')
            lines.append(line)

    # write file
    _writelines(path + '/' + file, lines)


### utility functions

def _split(str, sep):
    n = find(str, sep)
    if n >= 0:
        return str[:n], str[n + len(sep):]
    else:
        return str, ''


def _merge(*parts):
    return ''.join(parts)


def _writelines(file, lines):
    """ Writes text file
    """
    with open(file, 'w') as f:
        f.writelines(lines)

