
import sys

from seisflows.tools import array, unix
from seisflows.tools.tools import exists, findpath
from seisflows.tools.seismic import getpar, setpar


### input file writers

def write_sources(coords, path='.', ws=1., suffix=''):
    """ Writes source information to text file
    """
    sx, sy, sz = coords

    filename = findpath('seisflows.plugins') + '/' + 'solver/specfem2d/SOURCE'
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


def smooth_legacy(input_path='', output_path='', parameters=[], span=0.):
        solver = sys.modules['seisflows_solver']
        PATH = sys.modules['seisflows_paths']

        if not exists(input_path):
            raise Exception

        if not exists(output_path):
            unix.mkdir(output_path)

        if solver.mesh_properties.nproc!=1:
            raise NotImplementedError

        # intialize arrays
        kernels = {}
        for key in parameters or solver.parameters:
            kernels[key] = []

        coords = {}
        for key in ['x', 'z']:
            coords[key] = []

        # read kernels
        for key in parameters or solver.parameters:
            kernels[key] += solver.io.read_slice(input_path, key+'_kernel', 0)

        if not span:
            return kernels

        # read coordinates
        for key in ['x', 'z']:
            coords[key] += solver.io.read_slice(PATH.MODEL_INIT, key, 0)

        mesh = array.stack(coords['x'][0], coords['z'][0])

        # apply smoother
        for key in parameters or solver.parameters:
            kernels[key] = [array.meshsmooth(kernels[key][0], mesh, span)]

        # write smooth kernels
        for key in parameters or solver.parameters:
            solver.io.write_slice(kernels[key][0], output_path, key+'_kernel', 0)

