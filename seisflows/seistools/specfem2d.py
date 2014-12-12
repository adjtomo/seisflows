
import glob as _glob
import string as _string
import numpy as _np

from seisflows.tools import unix
from seisflows.tools.code import Struct
from seisflows.tools.config import findpath
from seisflows.seistools.segy import segyreader, segywriter


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
        i +=1

    # keep track of file names
    h.files = []
    for file in files:
        file = unix.basename(file)
        h.files.append(file)

    return s,h


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
        elif channel == ['p']:
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
            elif channel == ['p']:
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

def readsu(channel=[],prefix='SEM',suffix='.su'):
    """ Reads Seismic Unix file
    """
    if suffix=='':
        suffix = '.su'

    if channel in ['x']:
        file = '%s/Ux_file_single%s' % (prefix,suffix)
    elif channel in ['y']:
        file = '%s/Uy_file_single%s' % (prefix,suffix)
    elif channel in ['z']:
        file = '%s/Uz_file_single%s' % (prefix,suffix)
    elif channel in ['p']:
        file = '%s/Up_file_single%s' % (prefix,suffix)
    else:
        Exception

    # read data from file
    d,h = segyreader.readsu(file)
    return d,h


def writesu(d,h,channel=[],prefix='SEM',suffix='.su.adj'):
    """ Writes Seismic Unix file
    """
    if suffix=='':
        suffix = '.su.adj'

    if channel in ['x']:
        file = '%s/Ux_file_single%s' % (prefix,suffix)
    elif channel in ['y']:
        file = '%s/Uy_file_single%s' % (prefix,suffix)
    elif channel in ['z']:
        file = '%s/Uz_file_single%s' % (prefix,suffix)
    elif channel in ['p']:
        file = '%s/Up_file_single%s' % (prefix,suffix)
    else:
        Exception

    # write data to file
    segywriter.writesu(file,d,h)


### input file writers

def write_sources(par,hdr,path='.',suffix=''):
    """ Writes source information to text file
    """
    file = findpath('seistools')+'/'+'specfem2d/SOURCE'
    with open(file,'r') as f:
        lines = f.readlines()

    file = path+'/'+'DATA/SOURCE'+suffix
    _writelines(file,lines)

    # adjust source coordinates
    setpar( 'xs', hdr.sx[0], file )
    setpar( 'zs', hdr.sz[0], file )
    setpar( 'ts', hdr.ts, file )

    # adjust source amplitude
    try:
        fs = float(getpar('factor',file))
        setpar('factor',str(fs*hdr.fs),file)
    except:
        pass

    # adjust source wavelet
    if 1:
        # Ricker wavelet
        setpar( 'time_function_type', 1, file )
    elif 0:
        # first derivative of Gaussian
        setpar( 'time_function_type', 2, file )
    elif 0:
        # Gaussian
        setpar( 'time_function_type', 3, file )
    elif 0:
        # Dirac
        setpar( 'time_function_type', 4, file )
    elif 0:
        # Heaviside
        setpar( 'time_function_type', 5, file )

    setpar( 'f0', par['F0'], file )


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


def write_parameters(par,version='git-devel'):
    """ Writes parameters to text file
    """
    # read template
    file = findpath('seistools')+'/'+'specfem2d/par-'+version
    with open(file,'r') as f:
        lines = f.readlines()
    lines[-1] = ' '.join(['1',str(par.NX),'1',str(par.NZ),'1'])

    # write parameter file
    file = 'DATA/Par_file'
    _writelines(file,lines)
    setpar( 'xmin',      str(par.XMIN) )
    setpar( 'xmax',      str(par.XMAX) )
    setpar( 'nx',        str(par.NX) )
    setpar( 'nt',        str(par.NT) )
    setpar( 'deltat',    str(par.DT) )
    setpar( 'nsources',  str(1)      )

    # write interfaces file
    file = 'DATA/interfaces.dat'
    lines = []
    lines.extend('2\n')
    lines.extend('2\n')
    lines.extend('%f %f\n' % (par.XMIN,par.ZMIN))
    lines.extend('%f %f\n' % (par.XMAX,par.ZMIN))
    lines.extend('2\n')
    lines.extend('%f %f\n' % (par.XMIN,par.ZMAX))
    lines.extend('%f %f\n' % (par.XMAX,par.ZMAX))
    lines.extend(str(par.NZ))
    _writelines(file,lines)


def getpar(key,file='DATA/Par_file',sep='='):
    """ Reads parameter from parfile
    """
    with open(file,'r') as f:

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
    """ Writes parameter to parfile
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


###

def interpolate():
    parts = self.load('DATA/model_velocity.dat_output')

    # interpolate material parameters
    x = np.array(parts['x'][:]).T
    z = np.array(parts['z'][:]).T
    meshcoords = np.column_stack((x,z))
    x = model['x'].flatten()
    z = model['z'].flatten()
    gridcoords = np.column_stack((x,z))
    for key in ['rho','vp','vs']:
        v = model[key].flatten()
        parts[key] = [griddata(gridcoords,v,meshcoords,'linear')]
    self.save('DATA/model_velocity.dat_input',parts)


def mesher(self):

    seistools.specfem2d.setpar('SIMULATION_TYPE', '1')
    seistools.specfem2d.setpar('SAVE_FORWARD', '.false.')

    seistools.specfem2d.setpar('assign_external_model', '.false.')
    seistools.specfem2d.setpar('READ_EXTERNAL_SEP_FILE', '.false.')
    nt = seistools.specfem2d.getpar('nt')
    seistools.specfem2d.setpar('nt',str(100))

    with open('log.mesher','w') as f:
        subprocess.call(self.mesher_binary,stdout=f)
        subprocess.call(self.solver_binary,stdout=f)

    seistools.specfem2d.setpar('assign_external_model', '.true.')
    seistools.specfem2d.setpar('READ_EXTERNAL_SEP_FILE', '.true.')
    seistools.specfem2d.setpar('nt',str(nt))


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

