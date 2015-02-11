
import copy as _copy
import glob as _glob
import string as _string
import numpy as _np

from seisflows.tools import unix
from seisflows.tools.code import Struct

from seisflows.seistools.segy import segyreader


def ascii_specfem2d(**kwargs):
    """ Reads seismic traces from text files
    """
    files = glob(solver='2d',**kwargs)
    t = _np.loadtxt(files[0])[:, 0]
    h = Struct()
    h['t0'] = t[0]
    h['nr'] = len(files)
    h['ns'] = 1
    h['dt'] = _np.mean(_np.diff(t))
    h['nt'] = len(t)

    # read data
    s = _np.zeros((h['nt'], h['nr']))
    i = 0
    for file in files:
        s[:, i] = _np.loadtxt(file)[:, 1]
        i += 1

    # keep track of file names
    h.files = []
    for file in files:
        file = unix.basename(file)
        h.files.append(file)

    return s, h


def su_specfem2d(channel=None, prefix='SEM', suffix='.su'):
    """ Reads Seismic Unix file
    """
    if suffix == '':
        suffix = '.su'

    if channel in ['x']:
        file = '%s/Ux_file_single%s' % (prefix, suffix)
    elif channel in ['y']:
        file = '%s/Uy_file_single%s' % (prefix, suffix)
    elif channel in ['z']:
        file = '%s/Uz_file_single%s' % (prefix, suffix)
    elif channel in ['p']:
        file = '%s/Up_file_single%s' % (prefix, suffix)
    else:
        raise ValueError('CHANNEL must be one of the following: x y z p')

    # read data from file
    d, h = segyreader.readsu(file)
    return d, h


def ascii_specfem3d(**kwargs):
    """ Reads seismic traces from text files
    """
    files = glob(solver='3d',**kwargs)
    t = _np.loadtxt(files[0])[:, 0]
    h = Struct()
    h['t0'] = t[0]
    h['nr'] = len(files)
    h['ns'] = 1
    h['dt'] = _np.mean(_np.diff(t))
    h['nt'] = len(t)

    # read data
    s = _np.zeros((h['nt'], h['nr']))
    i = 0
    for file in files:
        s[:, i] = _np.loadtxt(file)[:, 1]
        i += 1

    # keep track of file names
    h.files = []
    for file in files:
        file = unix.basename(file)
        h.files.append(file)

    return s, h


def su_specfem3d(channel=None, prefix='SEM', suffix='', verbose=False):
    """ Reads Seismic Unix file
    """
    if channel in ['x']:
        wildcard = '%s/*_dx_SU%s' % (prefix, suffix)
    elif channel in ['y']:
        wildcard = '%s/*_dy_SU%s' % (prefix, suffix)
    elif channel in ['z']:
        wildcard = '%s/*_dz_SU%s' % (prefix, suffix)
    elif channel in ['p']:
        wildcard = '%s/*_dp_SU%s' % (prefix, suffix)
    else:
        raise ValueError('CHANNEL must be one of the following: x y z p')

    files = _glob.glob(wildcard)
    files = sorted(files, key=lambda x: int(unix.basename(x).split('_')[0]))

    file = files.pop(0)
    d, h = segyreader.readsu(file)

    if verbose:
        print file
        print 'number of traces:', d.shape[1]
        print 'min, max:', d.min(), d.max()
        print ''

    rx = _list(h.rx)
    ry = _list(h.ry)
    rz = _list(h.rz)
    sx = _list(h.sx)
    sy = _list(h.sy)
    sz = _list(h.sz)

    nn = [h.nr]
    nr = h.nr

    for file in files:
        d_, h_ = segyreader.readsu(file)

        # combine arrays
        d = _np.column_stack((d, d_))

        if verbose:
            print file
            print 'number of traces:', d_.shape[1]
            print 'min, max:', d_.min(), d_.max()
            print ''

        # combine headers
        rx.extend(h_.rx)
        ry.extend(h_.ry)
        rz.extend(h_.rz)
        sx.extend(h_.sx)
        sy.extend(h_.sy)
        sz.extend(h_.sz)
        nn.append(h_.nr)
        nr = nr + h_.nr

    h.rx = _np.array(rx)
    h.ry = _np.array(ry)
    h.rz = _np.array(rz)
    h.sx = _np.array(sx)
    h.sy = _np.array(sy)
    h.sz = _np.array(sz)

    h.nn = nn
    h.nr = nr

    return d, h


def ascii_specfem3d_globe(**kwargs):
    """ Reads seismic traces from text files
    """
    files = glob(solver='3d_globe', suffix='sem.ascii', **kwargs)
    t = _np.loadtxt(files[0])[:, 0]
    h = Struct()
    h['t0'] = t[0]
    h['nr'] = len(files)
    h['ns'] = 1
    h['dt'] = _np.mean(_np.diff(t))
    h['nt'] = len(t)

    # read data
    s = _np.zeros((h['nt'], h['nr']))
    i = 0
    for file in files:
        s[:, i] = _np.loadtxt(file)[:, 1]
        i += 1

    # keep track of file names
    h.files = []
    for file in files:
        file = unix.basename(file)
        h.files.append(file)

    return s, h


### utility functions

def glob(files=None, channel=None, prefix='SEM', suffix='semd', solver='3d'):
    """ Looks for seismic traces in current directory
    """
    if files:
        return files

    if solver=='2d':
        if channel in ['x']:
            wildcard = '%s/*.?XX.%s' % (prefix, suffix)
        elif channel in ['y']:
            wildcard = '%s/*.?XY.%s' % (prefix, suffix)
        elif channel in ['z']:
            wildcard = '%s/*.?XZ.%s' % (prefix, suffix)
        else:
            wildcard = '%s/*.?X?.%s' % (prefix, suffix)
    elif solver=='3d':
        if channel in ['x']:
            wildcard = '%s/*[xX]*.%s' % (prefix, suffix)
        elif channel in ['y']:
            wildcard = '%s/*[yY]*.%s' % (prefix, suffix)
        elif channel in ['z']:
            wildcard = '%s/*[zZ]*.%s' % (prefix, suffix)
        else:
            raise ValueError('CHANNEL must be one of the following: x y z p')
    elif solver=='3d_globe':
        if channel in ['e']:
            wildcard = '%s/*.?XE.%s' % (prefix, suffix)
        elif channel in ['n']:
            wildcard = '%s/*.?XN.%s' % (prefix, suffix)
        elif channel in ['z']:
            wildcard = '%s/*.?XZ.%s' % (prefix, suffix)
        else:
            wildcard = '%s/*.?X?.%s' % (prefix, suffix)


    files = _glob.glob(wildcard)
    if not files:
        raise Exception

    files.sort()
    return files


def _list(array):
    array2 = (_copy.copy(array))
    return list(array2)

