
import string as _string
import numpy as _np

from seisflows.tools import unix
from seisflows.tools.code import Struct
from seisflows.tools.config import findpath


### input file writers

def write_sources(par, hdr, path='.', suffix=''):
    """ Writes source information to text file
    """
    file = findpath('seistools') + '/' + 'specfem2d/SOURCE'
    with open(file, 'r') as f:
        lines = f.readlines()

    file = path + '/' + 'DATA/SOURCE' + suffix
    _writelines(file, lines)

    # adjust source coordinates
    setpar('xs', hdr.sx[0], file)
    setpar('zs', hdr.sz[0], file)
    setpar('ts', hdr.ts, file)

    # adjust source amplitude
    try:
        fs = float(getpar('factor', file))
        setpar('factor', str(fs*hdr.fs), file)
    except:
        pass

    # adjust source wavelet
    if 1:
        # Ricker wavelet
        setpar('time_function_type', 1, file)
    elif 0:
        # first derivative of Gaussian
        setpar('time_function_type', 2, file)
    elif 0:
        # Gaussian
        setpar('time_function_type', 3, file)
    elif 0:
        # Dirac
        setpar('time_function_type', 4, file)
    elif 0:
        # Heaviside
        setpar('time_function_type', 5, file)

    setpar('f0', par['F0'], file)


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


def write_parameters(par, version='git-devel'):
    """ Writes parameters to text file
    """
    # read template
    file = findpath('seistools') + '/' + 'specfem2d/par-' + version
    with open(file, 'r') as f:
        lines = f.readlines()
    lines[-1] = ' '.join(['1', str(par.NX), '1', str(par.NZ), '1'])

    # write parameter file
    file = 'DATA/Par_file'
    _writelines(file, lines)
    setpar('xmin', str(par.XMIN))
    setpar('xmax', str(par.XMAX))
    setpar('nx', str(par.NX))
    setpar('nt', str(par.NT))
    setpar('deltat', str(par.DT))
    setpar('nsources', str(1))

    # write interfaces file
    file = 'DATA/interfaces.dat'
    lines = []
    lines.extend('2\n')
    lines.extend('2\n')
    lines.extend('%f %f\n'%(par.XMIN, par.ZMIN))
    lines.extend('%f %f\n'%(par.XMAX, par.ZMIN))
    lines.extend('2\n')
    lines.extend('%f %f\n'%(par.XMIN, par.ZMAX))
    lines.extend('%f %f\n'%(par.XMAX, par.ZMAX))
    lines.extend(str(par.NZ))
    _writelines(file, lines)


def getpar(key, file='DATA/Par_file', sep='='):
    """ Reads parameter from parfile
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
    """ Writes parameter to parfile
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
    print rank
    return [names[ii] for ii in rank]


def _stack(a1, a2):
    if a1.size > 0:
        return _np.column_stack((a1, a2))
    else:
        return a2

