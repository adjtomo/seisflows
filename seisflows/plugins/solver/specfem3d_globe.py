
from seisflows.tools.tools import findpath
from seisflows.tools.seismic import setpar


def write_sources(PAR, h, path='.'):
    """ Writes source information to text file
    """
    filename = findpath('sesiflows.plugins') + '/' + 'specfem3d/SOURCE'
    with open(filename, 'r') as f:
        lines = f.readlines()

    filename = 'DATA/SOURCE'
    with open(filename, 'w') as f:
        f.writelines(lines)

    # adjust coordinates
    setpar('xs', h.sx[0], filename)
    setpar('zs', h.sz[0], filename)
    setpar('ts', h.ts, filename)

    # adjust wavelet
    setpar('f0', PAR['F0'], filename)


def write_receivers(h):
    """ Writes receiver information to text file
    """
    filename = 'DATA/STATIONS'
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

    with open(filename, 'w') as f:
        f.writelines(lines)


def write_parameters(par, version):
    """ Writes parameters to text file
    """
    raise NotImplementedError


