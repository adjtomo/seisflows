
from seisflows.tools.code import findpath
from seisflows.seistools.shared import getpar, setpar


### input file writers

def write_sources(coords, path='.', ws=1., suffix=''):
    """ Writes source information to text file
    """
    sx, sy, sz = coords

    file = findpath('seisflows.seistools') + '/' + 'specfem2d/SOURCE'
    with open(file, 'r') as f:
        lines = f.readlines()

    file = 'DATA/SOURCE' + suffix
    _writelines(file, lines)

    # adjust source coordinates
    setpar('xs', sx, file)
    setpar('zs', sy, file)
    #setpar('ts', ts[0], file)

    # adjust source amplitude
    try:
        fs = float(getpar('factor', file))
        fs *= ws
        setpar('factor', str(fs), file)
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

    #setpar('f0', par['F0'], file)


def write_receivers(coords, path='.'):
    """ Writes receiver information to text file
    """
    rx, ry, rz = coords
    nr = len(coords[0])

    filename = path +'/'+ 'DATA/STATIONS'

    lines = []
    for ir in range(nr):
        line = ''
        line += 'S%06d' % ir + ' '
        line += 'AA' + ' '
        line += '%11.5e' % rx[ir] + ' '
        line += '%11.5e' % ry[ir] + ' '
        line += '%3.1f' % 0. + ' '
        line += '%3.1f' % 0. + '\n'
        lines.extend(line)

    _writelines(filename, lines)


### utility functions

def _writelines(file, lines):
    """ Writes text file
    """
    with open(file, 'w') as f:
        f.writelines(lines)

