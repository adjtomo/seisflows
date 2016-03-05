
from seisflows.tools.code import findpath
from seisflows.seistools.shared import setpar



### input file writers

def write_sources(PAR, h, path='.'):
    """ Writes source information to text file
    """
    file = findpath('sesiflows.seistools') + '/' + 'specfem3d/SOURCE'
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


### utility functions

def _writelines(file, lines):
    """ Writes text file
    """
    with open(file, 'w') as f:
        f.writelines(lines)

