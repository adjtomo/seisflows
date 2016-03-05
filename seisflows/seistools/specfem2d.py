
from seisflows.tools.code import findpath
from seisflows.seistools.shared import getpar, setpar


### input file writers

def write_sources(par, hdr, path='.', suffix=''):
    """ Writes source information to text file
    """
    file = findpath('sesiflows.seistools') + '/' + 'specfem2d/SOURCE'
    with open(file, 'r') as f:
        lines = f.readlines()

    file = path + '/' + 'DATA/SOURCE' + suffix
    _writelines(file, lines)

    # adjust source coordinates
    setpar('xs', hdr.sx[0], file)
    setpar('zs', hdr.sy[0], file)
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
        line += '%11.5e' % h.ry[ir] + ' '
        line += '%3.1f' % 0. + ' '
        line += '%3.1f' % 0. + '\n'
        lines.extend(line)

    # write file
    _writelines(file, lines)


def write_parameters(par, version='git-devel'):
    """ Writes parameters to text file
    """
    # read template
    file = findpath('sesiflows.seistools') + '/' + 'specfem2d/par-' + version
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



### utility functions

def _writelines(file, lines):
    """ Writes text file
    """
    with open(file, 'w') as f:
        f.writelines(lines)

