import copy as _copy
import glob as _glob
import string as _string
import numpy as _np

from seisflows.tools import unix
from seisflows.tools.code import Struct
from seisflows.tools.config import findpath

import segy.reader as segyreader
import segy.writer as segywriter


### input file writers

def write_sources(PAR, h, path='.'):
    """ Writes source information to text file
    """
    file = findpath('seistools') + '/' + 'specfem3d/SOURCE'
    with open(file, 'r') as f:
        lines = f.readlines()

    file = 'DATA/SOURCE'
    _writelines(file, lines)

    # adjust coordinates
    setpar('xs', h.sx[0], file)
    setpar('zs', h.sz[0], file)
    setpar('ts', h.ts, file)

    # adjust wavelet
    setpar('f0', PAR['F0'], file)


def write_receivers(h):
    """ Writes receiver information to text file
    """
    file = 'DATA/STATIONS'
    lines = []

    # loop over receivers
    for ir in range(h.nr):
        line = ''
        line += 'S%06d' % ir + ' '
        line += 'AA' + ' '
        line += '%11.5e' % h.rx[ir] + ' '
        line += '%11.5e' % h.rz[ir] + ' '
        line += '%3.1f' % 0. + ' '
        line += '%3.1f' % 0. + '\n'
        lines.extend(line)

    # write file
    _writelines(file, lines)


def write_parameters(par, version):
    """ Writes parameters to text file
    """
    raise NotImplementedError


def getpar(key, file='DATA/Par_file', sep='='):
    """ Reads parameter from SPECFEM parfile
    """
    with open(file, 'r') as f:

        # read line by line
        for line in f:
            if _string.find(line, key) == 0:
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
            if _string.find(line, key) == 0:
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

def _writelines(file, lines):
    """ Writes text file
    """
    with open(file, 'w') as f:
        f.writelines(lines)


def _split(str, sep):
    from string import find

    n = find(str, sep)
    if n >= 0:
        return str[:n], str[n + len(sep):]
    else:
        return str, ''


def _merge(*parts):
    return ''.join(parts)


def _cmp(names):
    rank = []
    for name in names:
        ii = int(unix.basename(name).split('_')[0])
        rank.append(ii)
    return [names[ii] for ii in rank]


def _stack(a1, a2):
    if a1.size > 0:
        return _np.column_stack((a1, a2))
    else:
        return a2
