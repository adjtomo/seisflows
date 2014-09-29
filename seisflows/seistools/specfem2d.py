
import glob as _glob
import string as _string
import numpy as _np

from seisflows.tools import unix
from seisflows.tools.codetools import Struct
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
    if channel == 1:
      fmt = '%s/S%s.AA.%sX.%s' % (prefix,'%04d',char,suffix)
    elif channel == 2:
      fmt = '%s/S%s.AA.%sY.%s' % (prefix,'%04d',char,suffix)
    elif channel == 3:
      fmt = '%s/S%s.AA.%sZ.%s' % (prefix,'%04d',char,suffix)
    elif channel == 4:
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
      if channel == 1:
        label = ''.join([parts[2][:-1],'X'])
      elif channel == 2:
        label = ''.join([parts[2][:-1],'Y'])
      elif channel == 3:
        label = ''.join([parts[2][:-1],'Z'])
      elif channel == 4:
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
      if channel == 1:  wildcard = '%s/*.?XX.%s' % (prefix,suffix)
      elif channel == 2:  wildcard = '%s/*.?XY.%s' % (prefix,suffix)
      elif channel == 3:  wildcard = '%s/*.?XZ.%s' % (prefix,suffix)
      else: wildcard = '%s/*.?X?.%s' % (prefix,suffix)
    else:
      if channel == 1:  wildcard = '%s/*[xX]*.%s' % (prefix,suffix)
      elif channel == 2:  wildcard = '%s/*[yY]*.%s' % (prefix,suffix)
      elif channel == 3:  wildcard = '%s/*[zZ]*.%s' % (prefix,suffix)
    files = _glob.glob(wildcard)
    if files:
      files.sort()
    else:
      raise Exception
    return files

  elif filetype == 'su':
    if channel == 1:  file = '%s/Ux_file_single.%s' % (prefix,suffix)
    elif channel == 2:  file = '%s/Uy_file_single.%s' % (prefix,suffix)
    elif channel == 3:  file = '%s/Uz_file_single.%s' % (prefix,suffix)
    return file


###  reading and writing seismograms Seismic Unix data

def readsu(channel=[],prefix='SEM',suffix='.bin'):
    """ Reads Seismic Unix file
    """
    if suffix=='':
      suffix = '.bin'

    if channel == 1:
      file = '%s/Ux_file_single%s' % (prefix,suffix)
    elif channel == 2:
      file = '%s/Uy_file_single%s' % (prefix,suffix)
    elif channel == 3:
      file = '%s/Uz_file_single%s' % (prefix,suffix)
    elif channel == 4:  
      file = '%s/Up_file_single%s' % (prefix,suffix)
    else:
      Exception

    # read data from file
    d,h = segyreader.readsu(file)
    return d,h


def writesu(d,h,channel=[],prefix='SEM',suffix='.bin.adj'):
    """ Writes Seismic Unix file
    """
    if suffix=='':
      suffix = '.bin.adj'

    if channel == 1:  
      file = '%s/Ux_file_single%s' % (prefix,suffix)
    elif channel == 2:  
      file = '%s/Uy_file_single%s' % (prefix,suffix)
    elif channel == 3:  
      file = '%s/Uz_file_single%s' % (prefix,suffix)
    elif channel == 4:
      file = '%s/Up_file_single%s' % (prefix,suffix)
    else:
      Exception

    # write data to file
    segywriter.writesu(file,d,h)


### input file writers

def writesrc(PAR,h,path='.'):
    """ Writes source information to text file
    """
    from seisflows.tools.configure import getpath

    file = getpath('seistools')+'/'+'specfem2d/SOURCE'
    with open(file,'r') as f:
      lines = f.readlines()

    file = 'DATA/SOURCE'
    _writelines(file,lines)

    # adjust source coordinates
    setpar( 'xs', h.sx[0], file )
    setpar( 'zs', h.sz[0], file )
    setpar( 'ts', h.ts, file )

    # adjust source amplitude
    try:
      fs = float(getpar('factor',file))
      setpar('factor',str(fs*h.fs),file)
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

    setpar( 'f0', PAR['F0'], file )


def writerec(nr,rx,rz):
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


def writepar(vars,version='git-devel'):
    """ Writes parameters to text file
    """
    from seisflows.tools.configure import getpath

    PAR = Struct(vars)

    # read template
    file = getpath('seistools')+'/'+'specfem2d/par-'+version
    with open(file,'r') as f:
      lines = f.readlines()
    lines[-1] = ' '.join(['1',str(PAR.NX),'1',str(PAR.NZ),'1'])

    # write parameter file
    file = 'DATA/Par_file'
    _writelines(file,lines)
    setpar( 'xmin',      str(PAR.XMIN) )
    setpar( 'xmax',      str(PAR.XMAX) )
    setpar( 'nx',        str(PAR.NX) )
    setpar( 'nt',        str(PAR.NT) )
    setpar( 'deltat',    str(PAR.DT) )
    setpar( 'nsources',  str(1)      )

    # write interfaces file
    file = 'DATA/interfaces.dat'
    lines = []
    lines.extend('2\n')
    lines.extend('2\n')
    lines.extend('%f %f\n' % (PAR.XMIN,PAR.ZMIN))
    lines.extend('%f %f\n' % (PAR.XMAX,PAR.ZMIN))
    lines.extend('2\n')
    lines.extend('%f %f\n' % (PAR.XMIN,PAR.ZMAX))
    lines.extend('%f %f\n' % (PAR.XMAX,PAR.ZMAX))
    lines.extend(str(PAR.NZ))
    _writelines(file,lines)


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
