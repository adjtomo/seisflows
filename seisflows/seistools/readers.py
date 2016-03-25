
import copy as _copy
import glob as _glob
import string as _string
import numpy as _np

from os.path import basename

from seisflows.tools.code import Struct

from seisflows.seistools.segy import segyreader


def ascii_specfem2d_obspy(**kwargs):
    """ Reads seismic traces from text files
    """
    from obspy.core.stream import Stream
    from obspy.core.trace import Trace

    filenames = glob(solver='specfem2d', **kwargs)

    t = _np.loadtxt(files[0])[:,0]
    nt = len(t)
    nr = len(filenames)

    d = Trace(data=np.zeros(nt, dtype='float32'))

    trace.stats.starttime = t[0]
    trace.stats.delta = _np.mean(_np.diff(t))
    trace.stats.nt = len(t)

    # read data
    stream = Stream(t)*nr

    for filename in filenames:
        stream.data = _np.loadtxt(filename)[:, 1]

    return stream


def ascii_specfem2d(**kwargs):
    """ Reads seismic traces from text files
    """
    files = glob(solver='specfem2d', **kwargs)
    t = _np.loadtxt(files[0])[:,0]
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
        file = basename(file)
        h.files.append(file)

    return s, h


def su_specfem2d_obspy(prefix='SEM', channel=None, suffix='.su'):
    from obspy import read

    if channel in ['x']:
        filename = '%s/Ux_file_single%s' % (prefix, suffix)
    elif channel in ['y']:
        filename = '%s/Uy_file_single%s' % (prefix, suffix)
    elif channel in ['z']:
        filename = '%s/Uz_file_single%s' % (prefix, suffix)
    elif channel in ['p']:
        filename = '%s/Up_file_single%s' % (prefix, suffix)
    else:
        raise ValueError('CHANNEL must be one of the following: x y z p')

    streamobj = read(filename, format='SU', byteorder='<')
    return streamobj


def su_specfem2d(prefix='SEM', channel=None, suffix='.su'):
    """ Reads Seismic Unix file
    """
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
    files = glob(solver='specfem3d',**kwargs)
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
        file = basename(file)
        h.files.append(file)

    return s, h


def su_specfem3d(prefix='SEM', channel=None, suffix='', verbose=False):
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
    files = sorted(files, key=lambda x: int(basename(x).split('_')[0]))

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

    # define proc number
    i_proc_buf = file.split('/')[-1]
    i_proc = int(i_proc_buf.split('_')[0])
    ip = [i_proc]

    for file in files:
        d_, h_ = segyreader.readsu(file)

        # define proc number  
        i_proc_buf = file.split('/')[-1]
        i_proc = int(i_proc_buf.split('_')[0])

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
        ip.append(i_proc)

    h.rx = _np.array(rx)
    h.ry = _np.array(ry)
    h.rz = _np.array(rz)
    h.sx = _np.array(sx)
    h.sy = _np.array(sy)
    h.sz = _np.array(sz)

    h.nn = nn
    h.nr = nr
    h.ip = ip

    return d, h


def su_specfem3d_obspy(prefix='SEM', channel=None, suffix='', byteorder='<', verbose=False):
    """ Reads Seismic Unix file
    """
    from obspy import read

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

    filenamess = _glob.glob(wildcard)

    sort_by = lambda x: int(basename(x).split('_')[0])
    filenames = sorted(filenamess, key=sort_by)

    streamobj = read(filenames.pop(), format='SU', byteorder='<')
    for filename in filenames:
        streamobj += read(filename, format='SU', byteorder='<')

    return streamobj


def ascii_specfem3d_globe(**kwargs):
    """ Reads seismic traces from text files
    """
    files = glob(solver='specfem3d_globe', suffix='sem.ascii', **kwargs)
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
        file = basename(file)
        h.files.append(file)

    return s, h


### utility functions

def glob(files=None, prefix='SEM', channel=None, suffix='semd', solver='specfem3d'):
    """ Looks for seismic traces in current directory
    """
    if files:
        return files

    if solver=='specfem2d':
        if channel in ['x']:
            wildcard = '%s/*.?XX.%s' % (prefix, suffix)
        elif channel in ['y']:
            wildcard = '%s/*.?XY.%s' % (prefix, suffix)
        elif channel in ['z']:
            wildcard = '%s/*.?XZ.%s' % (prefix, suffix)
        else:
            wildcard = '%s/*.?X?.%s' % (prefix, suffix)
    elif solver=='specfem3d':
        if channel in ['x']:
            wildcard = '%s/*[xX]*.%s' % (prefix, suffix)
        elif channel in ['y']:
            wildcard = '%s/*[yY]*.%s' % (prefix, suffix)
        elif channel in ['z']:
            wildcard = '%s/*[zZ]*.%s' % (prefix, suffix)
        else:
            raise ValueError('CHANNEL must be one of the following: x y z p')
    elif solver=='specfem3d_globe':
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

