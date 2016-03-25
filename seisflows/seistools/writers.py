
import string as _string
import numpy as _np

from seisflows.tools import unix

from seisflows.seistools.shared import SeisStruct
from seisflows.seistools.segy import segywriter


def ascii_specfem2d(f, h, prefix='SEM', channel=None, suffix='adj', char='FX', opt=''):
    """ Writes seismic traces to text files
    """

    files = []

    if opt == 'legacy':
        if channel in ['x']:
            fmt = '%s/S%s.AA.%sX.%s' % (prefix, '%04d', char, suffix)
        elif channel in ['y']:
            fmt = '%s/S%s.AA.%sY.%s' % (prefix, '%04d', char, suffix)
        elif channel in ['z']:
            fmt = '%s/S%s.AA.%sZ.%s' % (prefix, '%04d', char, suffix)
        elif channel == ['p']:
            fmt = '%s/S%s.AA.%sP.%s' % (prefix, '%04d', char, suffix)
        else:
            raise ValueError('CHANNEL must be one of the following: x y z p')
        for i in range(h.nr):
            files.append(fmt % (i+1))

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin + h['nt'])
        t = _np.arange(imin, imax)*h['dt']

        for i in range(h.nr):
            w = f[:, i]
            _np.savetxt(files[i], _np.column_stack((t, w)), '%11.4e')

    else:
        for file in h.files:
            parts = _string.split(file, '.')
            if channel in ['x']:
                label = ''.join([parts[2][:-1], 'X'])
            elif channel in ['y']:
                label = ''.join([parts[2][:-1], 'Y'])
            elif channel in ['z']:
                label = ''.join([parts[2][:-1], 'Z'])
            elif channel == ['p']:
                label = ''.join([parts[2][:-1], 'P'])
            else:
                raise ValueError('CHANNEL must be one of the following: x y z p')

            parts[-2] = label
            parts[-1] = 'adj'

            files.append(prefix + '/' + '.'.join(parts))

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin + h['nt'])
        t = _np.arange(imin, imax)*h['dt']

        for i, file in enumerate(files):
            w = f[:, i]
            _np.savetxt(file, _np.column_stack((t, w)), '%11.4e')


def su_specfem2d_obspy(d, prefix='SEM', channel=None, suffix='.su.adj'):
    """ Writes Seismic Unix file
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

    max_delta = 0.065535
    dummy_delta = max_delta

    if d[0].stats.delta > max_delta:
        for t in d:
            t.stats.delta = dummy_delta

    # write data to file
    d.write(file, format='SU')



def su_specfem2d(d, h, prefix='SEM', channel=None, suffix='.su.adj'):
    """ Writes Seismic Unix file
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

    # write data to file
    segywriter.writesu(file, d, h)


def ascii_specfem3d(f, h, prefix='SEM', channel=None, suffix='adj', char='FX', opt=''):
    """ Writes seismic traces to text files
    """

    files = []

    if opt == 'legacy':
        if channel in ['x']:
            fmt = '%s/S%s.AA.%sX.%s' % (prefix, '%04d', char, suffix)
        elif channel in ['y']:
            fmt = '%s/S%s.AA.%sY.%s' % (prefix, '%04d', char, suffix)
        elif channel in ['z']:
            fmt = '%s/S%s.AA.%sZ.%s' % (prefix, '%04d', char, suffix)
        elif channel in ['p']:
            fmt = '%s/S%s.AA.%sP.%s' % (prefix, '%04d', char, suffix)
        else:
            raise ValueError('CHANNEL must be one of the following: x y z p')
        for i in range(h.nr):
            files.append(fmt % (i + 1))

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin + h['nt'])
        t = _np.arange(imin, imax)*h['dt']

        for i in range(h.nr):
            w = f[:, i]
            _np.savetxt(files[i], _np.column_stack((t, w)), '%11.4e')

    else:
        for file in h.files:
            parts = _string.split(file, '.')
            if channel in ['x']:
                label = ''.join([parts[2][:-1], 'X'])
            elif channel in ['y']:
                label = ''.join([parts[2][:-1], 'Y'])
            elif channel in ['z']:
                label = ''.join([parts[2][:-1], 'Z'])
            elif channel in ['p']:
                label = ''.join([parts[2][:-1], 'P'])

            parts[-2] = label
            parts[-1] = 'adj'

            files.append(prefix + '.'.join(parts))

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin + h['nt'])
        t = _np.arange(imin, imax)*h['dt']

        for i, file in enumerate(files):
            w = f[:, i]
            _np.savetxt(file, _np.column_stack((t, w)), '%11.4e')


def su_specfem3d(d, h, prefix='SEM', channel=None, suffix='.adj', verbose=False):
    nproc = len(h.nn)

    if channel in ['x']:
        wildcard = '%s/%d_dx_SU%s'
    elif channel in ['y']:
        wildcard = '%s/%d_dy_SU%s'
    elif channel in ['z']:
        wildcard = '%s/%d_dz_SU%s'
    elif channel in ['p']:
        wildcard = '%s/%d_dp_SU%s'
    else:
        raise ValueError('CHANNEL must be one of the following: x y z p')


    imax = 0

    for iproc in range(nproc):

        file = wildcard % (prefix, h.ip[iproc], suffix)
        imin = imax
        imax = imax + h.nn[iproc]

        d_ = d[:, imin:imax]
        h_ = SeisStruct(nr=h.nr, nt=h.nt, dt=h.dt, ts=h.ts, nrec=h.nrec,
                        nsrc=h.nsrc)
        h_.rx = h.rx[imin:imax]
        h_.ry = h.ry[imin:imax]
        h_.rz = h.rz[imin:imax]
        h_.sx = h.sx[imin:imax]
        h_.sy = h.sy[imin:imax]
        h_.sz = h.sz[imin:imax]

        h_.nr = imax - imin

        if verbose:
            print file
            print (imin, imax)
            print ''

        segywriter.writesu(file, d_, h_)


def ascii_specfem3d_globe(f, h, prefix='SEM', channel=None, suffix='adj', char='FX', opt=''):
    """ Writes seismic traces to text files
    """

    files = []

    if opt == 'legacy':
        if channel in ['x']:
            fmt = '%s/S%s.AA.%sE.%s' % (prefix, '%04d', char, suffix)
        elif channel in ['y']:
            fmt = '%s/S%s.AA.%sN.%s' % (prefix, '%04d', char, suffix)
        elif channel in ['z']:
            fmt = '%s/S%s.AA.%sZ.%s' % (prefix, '%04d', char, suffix)
        elif channel in ['p']:
            fmt = '%s/S%s.AA.%sP.%s' % (prefix, '%04d', char, suffix)
        else:
            raise ValueError('CHANNEL must be one of the following: x y z p')
        for i in range(h.nr):
            files.append(fmt % (i + 1))

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin + h['nt'])
        t = _np.arange(imin, imax)*h['dt']

        for i in range(h.nr):
            w = f[:, i]
            _np.savetxt(files[i], _np.column_stack((t, w)), '%11.4e')

    else:
        for file in h.files:
            parts = _string.split(file, '.')[:-1]
            if channel in ['x']:
                label = ''.join([parts[2][:-1], 'E'])
            elif channel in ['y']:
                label = ''.join([parts[2][:-1], 'N'])
            elif channel in ['z']:
                label = ''.join([parts[2][:-1], 'Z'])
            elif channel in ['p']:
                label = ''.join([parts[2][:-1], 'P'])
            else:
                raise ValueError('CHANNEL must be one of the following: x y z p')

            parts[-2] = label
            parts[-1] = 'adj'

            files.append(prefix + '/' + '.'.join(parts))

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin + h['nt'])
        t = _np.arange(imin, imax)*h['dt']

        for i, file in enumerate(files):
            w = f[:, i]
            _np.savetxt(file, _np.column_stack((t, w)), '%11.4e')

