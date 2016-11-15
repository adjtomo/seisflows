
from seisflows.tools.code import findpath
from seisflows.tools.shared import getpar, setpar


### input file writers

def write_sources(coords, path='.', ws=1., suffix=''):
    """ Writes source information to text file
    """
    sx, sy, sz = coords

    filename = findpath('seisflows.plugins') + '/' + 'specfem2d/SOURCE'
    with open(filename, 'r') as f:
        lines = f.readlines()

    filename = 'DATA/SOURCE' + suffix
    with open(filename, 'w') as f:
        f.writelines(lines)

    # adjust source coordinates
    setpar('xs', sx, filename)
    setpar('zs', sy, filename)
    #setpar('ts', ts[0], filename)

    # adjust source amplitude
    try:
        fs = float(getpar('factor', filename))
        fs *= ws
        setpar('factor', str(fs), filename)
    except:
        pass

    # adjust source wavelet
    if 1:
        # Ricker wavelet
        setpar('time_function_type', 1, filename)
    elif 0:
        # first derivative of Gaussian
        setpar('time_function_type', 2, filename)
    elif 0:
        # Gaussian
        setpar('time_function_type', 3, filename)
    elif 0:
        # Dirac
        setpar('time_function_type', 4, filename)
    elif 0:
        # Heaviside
        setpar('time_function_type', 5, filename)

    #setpar('f0', par['F0'], filename)


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

    with open(filename, 'w') as f:
        f.writelines(lines)


