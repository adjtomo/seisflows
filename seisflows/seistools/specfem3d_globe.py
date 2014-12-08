
import copy as _copy
import glob as _glob
import string as _string
import numpy as _np

from seisflows.tools import unix
from seisflows.tools.codetools import Struct
from seisflows.tools.configtools import findpath

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
            fmt = '%s/S%s.AA.%sE.%s' % (prefix,'%04d',char,suffix)
        elif channel in ['y']:
            fmt = '%s/S%s.AA.%sN.%s' % (prefix,'%04d',char,suffix)
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
            parts = _string.split(file,'.')[:-1]
            if channel in ['x']:
                label = ''.join([parts[2][:-1],'E'])
            elif channel in ['y']:
                label = ''.join([parts[2][:-1],'N'])
            elif channel in ['z']:
                label = ''.join([parts[2][:-1],'Z'])
            elif channel in ['p']:
                label = ''.join([parts[2][:-1],'P'])

            parts[-2] = label
            parts[-1] = 'adj'

            files.append( prefix+'/'+'.'.join(parts) )

        # write data to files
        imin = int(_np.floor(h['t0']/h['dt']))
        imax = int(imin+h['nt'])
        t = _np.arange(imin,imax)*h['dt']

        for i,file in enumerate(files):
            w = f[:,i]
            _np.savetxt(file,_np.column_stack((t,w)),'%11.4e')


def glob(files=[],filetype='ascii',channel=[],prefix='SEM',suffix='sem.ascii'):
    """ Checks for seismic traces in current directory
    """
    if files:
        return files

    elif filetype == 'ascii':
        if 'sem' in suffix:
            if channel in ['e']:  wildcard = '%s/*.?XE.%s' % (prefix,suffix)
            elif channel in ['n']:  wildcard = '%s/*.?XN.%s' % (prefix,suffix)
            elif channel in ['z']:  wildcard = '%s/*.?XZ.%s' % (prefix,suffix)
            else: wildcard = '%s/*.?X?.%s' % (prefix,suffix)
        files = _glob.glob(wildcard)
        if files:
            files.sort()
        else:
            raise Exception
        return files


### input file writers

def write_sources(PAR,h,path='.'):
    """ Writes source information to text file
    """
    file = findpath('seistools')+'/'+'specfem3d/SOURCE'
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



def write_receivers(h):
    """ Writes receiver information to text file
    """
    file = 'DATA/STATIONS'
    lines = []

    # loop over receivers
    for ir in range(h.nr):
        line = ''
        line += 'S%06d'  % ir       + ' '
        line += 'AA'                + ' '
        line += '%11.5e' % h.rx[ir] + ' '
        line += '%11.5e' % h.rz[ir] + ' '
        line += '%3.1f'  % 0.       + ' '
        line += '%3.1f'  % 0.       + '\n'
        lines.extend(line)

    # write file
    _writelines(file,lines)


def write_parameters(par,version):
    """ Writes parameters to text file
    """
    raise NotImplementedError


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
    return [names[ii] for ii in rank]


def _stack(a1,a2):
    if a1.size > 0:
        return _np.column_stack((a1,a2))
    else:
        return a2
