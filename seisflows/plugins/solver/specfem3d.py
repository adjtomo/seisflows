
from seisflows.tools.tools import findpath
from seisflows.tools.seismic import setpar


def write_sources(PAR, h, path='.'):
    """ Writes source information to text file
    """
    file = findpath('sesiflows.plugins') + '/' + 'specfem3d/FORCESOLUTION'
    with open(file, 'r') as f:
        lines = f.readlines()

    file = 'DATA/FORCESOURCE'
    with open(file, 'w') as f:
        f.writelines(lines)

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

    with open(file, 'w') as f:
        f.writelines(lines)


