
import copy as _copy
import glob as _glob
import string as _string
import numpy as _np

from seisflows.tools import unix
from seisflows.tools.codetools import Struct

import segy.reader as segyreader
import segy.writer as segywriter


### reading and writing ASCII data

def read(**kwargs):
    """ Reads seismic traces from text files
    """
    files = glob(**kwargs)
    t = _np.loadtxt(files[0])[:,0]
    h = Struct()
    h['t0'] = t[0]
    h['nr'] = len(files)
    h['ns'] = 1
    h['dt'] = _np.mean(_np.diff(t))
    h['nt'] = len(t)

    # read data
    s = _np.zeros((h['nt'],h['nr']))
    i = 0
    for file in files:
        s[:,i] = _np.loadtxt(file)[:,1]
        i = i+1

    # keep track of file names
    h.files = []
    for file in files:
        file = unix.basename(file)
        h.files.append(file)

    return (s,h)


def write(f,h,channel,char='FX',prefix='SEM',suffix='adj',opt=''):
    """ Writes seismic traces to text files
    """

    files = []

    if opt=='legacy':
        if channel in ['x']:
            fmt = '%s/S%s.AA.%sX.%s' % (prefix,'%04d',char,suffix)
        elif channel in ['y']:
            fmt = '%s/S%s.AA.%sY.%s' % (prefix,'%04d',char,suffix)
        elif channel in ['z']:
            fmt = '%s/S%s.AA.%sZ.%s' % (prefix,'%04d',char,suffix)
        elif channel in ['p']:
            fmt = '%s/S%s.AA.%sP.%s' % (prefix,'%04d',char,suffix)
        for i in range(h.nr):
            files.append( fmt % (i+1) )

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin+h['nt'])
        t = _np.arange(imin,imax)*h['dt']

        for i in range(h.nr):
            w = f[:,i]
            _np.savetxt(files[i],_np.column_stack((t,w)),'%11.4e')

    else:
        for file in h.files:
            parts = _string.split(file,'.')
            if channel in ['x']:
                label = ''.join([parts[2][:-1],'X'])
            elif channel in ['y']:
                label = ''.join([parts[2][:-1],'Y'])
            elif channel in ['z']:
                label = ''.join([parts[2][:-1],'Z'])
            elif channel in ['p']:
                label = ''.join([parts[2][:-1],'P'])

            parts[-2] = label
            parts[-1] = 'adj'

            files.append( prefix+'.'.join(parts) )

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin+h['nt'])
        t = _np.arange(imin,imax)*h['dt']

        for i,file in enumerate(files):
            w = f[:,i]
            _np.savetxt(file,_np.column_stack((t,w)),'%11.4e')


def glob(files=[],filetype='ascii',channel=[],prefix='SEM',suffix='semd'):
    """ Checks for seismic traces in current directory
    """
    if files:
        return files

    elif filetype == 'ascii':
        if 'sem' in suffix:
            if channel in ['x']:  wildcard = '%s/*.?XX.%s' % (prefix,suffix)
            elif channel in ['y']:  wildcard = '%s/*.?XY.%s' % (prefix,suffix)
            elif channel in ['z']:  wildcard = '%s/*.?XZ.%s' % (prefix,suffix)
            else: wildcard = '%s/*.?X?.%s' % (prefix,suffix)
        else:
            if channel in ['x']:  wildcard = '%s/*[xX]*.%s' % (prefix,suffix)
            elif channel in ['y']:  wildcard = '%s/*[yY]*.%s' % (prefix,suffix)
            elif channel in ['z']:  wildcard = '%s/*[zZ]*.%s' % (prefix,suffix)
        files = _glob.glob(wildcard)
        if files:
            files.sort()
        else:
            raise Exception
        return files

    elif filetype == 'su':
        if channel in ['x']:  file = '%s/Ux_file_single.%s' % (prefix,suffix)
        elif channel in ['y']:  file = '%s/Uy_file_single.%s' % (prefix,suffix)
        elif channel in ['z']:  file = '%s/Uz_file_single.%s' % (prefix,suffix)
        return file


###  reading and writing seismograms Seismic Unix data

def readsu(channel=[],prefix='SEM',suffix='',verbose=False):
    """ Reads Seismic Unix file
    """
    if channel in ['x']:
        wildcard = '%s/*_dx_SU%s' % (prefix,suffix)
    elif channel in ['y']:
        wildcard = '%s/*_dy_SU%s' % (prefix,suffix)
    elif channel in ['z']:
        wildcard = '%s/*_dz_SU%s' % (prefix,suffix)
    elif channel in ['p']:
        wildcard = '%s/*_dp_SU%s' % (prefix,suffix)
    else:
        Exception

    files = _glob.glob(wildcard)
    files = sorted(files,key=lambda x:int(unix.basename(x).split('_')[0]))

    file = files.pop(0)
    d,h = segyreader.readsu(file)

    if verbose:
        print file
        print 'number of traces:', d.shape[1]
        print 'min, max:', d.min(), d.max()
        print ''

    rx = list(_copy.copy(h.rx))
    ry = list(_copy.copy(h.ry))
    rz = list(_copy.copy(h.rz))
    sx = list(_copy.copy(h.sx))
    sy = list(_copy.copy(h.sy))
    sz = list(_copy.copy(h.sz))

    nn = [h.nr]
    nr = h.nr

    for file in files:
        d_,h_ = segyreader.readsu(file)

        # combine arrays
        d = _np.column_stack((d,d_))

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

    return d,h


def writesu(d,h,channel=[],prefix='SEM',suffix='.adj',verbose=False):

    nproc = len(h.nn)

    if suffix=='':
        suffix = '.adj'

    if channel in ['x']:
        wildcard = '%s/%d_dx_SU%s'
    elif channel in ['y']:
        wildcard = '%s/%d_dy_SU%s'
    elif channel in ['z']:
        wildcard = '%s/%d_dz_SU%s'
    elif channel in ['p']:
        wildcard = '%s/%d_dp_SU%s'
    else:
        Exception

    imin=0
    imax=0

    for iproc in range(nproc):

        file = wildcard % (prefix,iproc,suffix)
        imin = imax
        imax = imax+h.nn[iproc]

        d_ = d[:,imin:imax]
        h_ = _copy.copy(h)
        h_.rx = h.rx[imin:imax]
        h_.ry = h.ry[imin:imax]
        h_.rz = h.rz[imin:imax]
        h_.sx = h.sx[imin:imax]
        h_.sy = h.sy[imin:imax]
        h_.sz = h.sz[imin:imax]

        h_.nr = imax-imin

        if verbose:
            print file
            print (imin,imax)
            print ''

        segywriter.writesu(file,d_,h_)


### input file writers

def write_sources(PAR,h,path='.'):
    """ Writes source information to text file
    """
    from seisflows.tools.configure import getpath

    file = getpath('seistools')+'/'+'specfem3d/SOURCE'
    with open(file,'r') as f:
        lines = f.readlines()

    file = 'DATA/SOURCE'
    _writelines(file,lines)

    # adjust coordinates
    setpar( 'xs', h.sx[0], file )
    setpar( 'zs', h.sz[0], file )
    setpar( 'ts', h.ts, file )

    # adjust wavelet
    setpar( 'f0', PAR['F0'], file )



def write_receivers(nr,rx,rz):
    """ Writes receiver information to text file
    """
    file = 'DATA/STATIONS'
    lines = []

    # loop over receivers
    for ir in range(nr):
        line = ''
        line = line + 'S%06d'  % ir     + ' '
        line = line + 'AA'              + ' '
        line = line + '%11.5e' % rx[ir] + ' '
        line = line + '%11.5e' % rz[ir] + ' '
        line = line + '%3.1f'  % 0.     + ' '
        line = line + '%3.1f'  % 0.     + '\n'
        lines.extend(line)

    # write file
    _writelines(file,lines)


def write_parameters(PAR):
    """ Writes parameters to text file
    """
    from seisflows.tools.configure import getpath

    PAR = Struct(PAR)

    # read template
    file = getpath('seistools')+'/'+'specfem3d/par-'+version
    with open(file,'r') as f:
        lines = f.readlines()

    # write parameter file

    # write interfaces file

    raise Exception


def getpar(key,file='DATA/Par_file',sep='='):
    """ Reads parameter from SPECFEM parfile
    """
    with open(file,'r') as f:
        lines = []

        # read line by line
        for line in f:
            if _string.find(line,key) == 0:
            # read key
                key,val = _split(line,sep)
                if not key:
                    continue
                # read val
                val,_ = _split(val,'#')
                return val.strip()

        raise Exception


def setpar(key,val,file='DATA/Par_file',path='.',sep='='):
    """ Writes parameter to SPECFEM parfile
    """

    val = str(val)

    # read line by line
    with open(path+'/'+file,'r') as f:
        lines = []
        for line in f:
            if _string.find(line,key) == 0:
                # read key
                key,_ = _split(line,sep)
                # read comment
                _,comment = _split(line,'#')
                n = len(line)-len(key)-len(val)-len(comment)-2
                # replace line
                if comment:
                    line = _merge(key,sep,val,' '*n,'#',comment)
                else:
                    line = _merge(key,sep,str(val),'\n')
            lines.append(line)

    # write file
    _writelines(path+'/'+file,lines)


### utility functions

def _writelines(file,lines):
    """ Writes text file
    """
    with open(file,'w') as f:
        f.writelines(lines)


def _split(str,sep):
    from string import find

    n = find(str,sep)
    if n >= 0:
        return str[:n], str[n+len(sep):]
    else:
        return str, ''


def _merge(*parts):
    return ''.join(parts)


def _cmp(names):
    rank = []
    for name in names:
        ii = int(unix.basename(name).split('_')[0])
        rank.append(ii)
    print rank
    return [names[ii] for ii in rank]


def _stack(a1,a2):
    if a1.size > 0:
        return _np.column_stack((a1,a2))
    else:
        return a2
