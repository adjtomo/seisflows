
import subprocess

import numpy as np

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configure import getclass, getpath, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class specfem2d(object):
  """ Python interface and base class for SPECFEM2D

    In the code, we distinguish between high-level and low-level interfaces.
    High-level methods deal with evaluation of the misfit function and its 
    derivatives and include 'evaluate_func', 'evaluate_grad', 'apply_hess'. 
    Low-level methods provide direct  accesss to the mesher, forward solver, 
    adjoint solver, and other SPECFEM2D components.
  """


  def __init__(self):
    """ Class constructor
    """
    # check user supplied parameters
    if 'XMIN' not in PAR or 'XMAX' not in PAR:
        raise Exception

    if 'ZMIN' not in PAR or 'ZMAX' not in PAR:
        raise Exception

    if 'NX' not in PAR:
        raise Exception

    if 'NZ' not in PAR:
        raise Exception

    if 'NT' not in PAR:
        raise Exception

    if 'DT' not in PAR:
        raise Exception

    if 'F0' not in PAR:
        raise Exception

    if 'WAVELET' not in PAR:
        setattr(PAR,'WAVELET','ricker')

    if 'PREPROCESS' not in PAR:
        setattr(PAR,'PREPROCESS','default')

    if 'MODEL' not in PAR:
        setattr(PAR,'MODEL','acoustic')

    # check user supplied paths
    if not exists(PATH.MODEL_INIT):
        raise Exception

    if not exists(PATH.DATA): 
        assert exists(PATH.MODEL_TRUE)
        assert exists(PATH.SOLVER_FILES)

    if not exists(PATH.SOLVER_BINARIES):
        raise Exception

    #if not exists(PATH.MESH):
    #    raise Exception

    # load preprocessing tools
    self.preprocess = getclass('preprocess',PAR.PREPROCESS)(
      reader=seistools.specfem2d.readsu,
      writer=seistools.specfem2d.writesu)

    # model parameters expected by solver
    model_parameters = []
    model_parameters += ['rho']
    model_parameters += ['vp']
    model_parameters += ['vs']
    self.model_parameters = model_parameters

    # model paramters included in inversion
    inversion_parameters = []
    if PAR.MODEL == 'acoustic':
      inversion_parameters += ['vs']
      self.inversion_parameters = inversion_parameters

    elif PAR.MODEL == 'elastic':
      inversion_parameters += ['rhop']
      inversion_parameters += ['vp']
      inversion_parameters += ['vs']
      self.inversion_parameters = inversion_parameters


  def setup(self,model_type='gll'):
      """ Prepares directories in which to run solver
      """
      unix.rm(self.getpath())
      unix.mkdir(self.getpath())
      unix.cd(self.getpath())

      # create subdirectories
      unix.mkdir('bin')
      unix.mkdir('DATA')
      unix.mkdir('OUTPUT_FILES')

      unix.mkdir('traces/obs')
      unix.mkdir('traces/syn')
      unix.mkdir('traces/adj')

      # copy binaries
      src = glob(PATH.SOLVER_BINARIES+'/'+'*')
      dst = 'bin/'
      unix.cp(src,dst)


      if PATH.DATA:
          # copy user supplied data
          src = glob(PATH.DATA+'/'+self.getshot()+'/'+'*')
          dst = 'traces/obs/'
          unix.cp(src,dst)

          # generate SPECFEM2D input files
          self.writepar()
          self.writesrc()
          self.writerec()

          # prepare starting model
          self.generate_mesh(
            model_path = PATH.MODEL_INIT,
            model_type = model_type,
            model_name = 'model_init')

      else:
         # copy user supplied SPECFEM2D input files
          src = glob(PATH.SOLVER_FILES+'/'+'*')
          dst = 'DATA/'
          unix.cp(src,dst)

          # generate data
          self.generate_data(
            model_path = PATH.MODEL_TRUE,
            model_type = model_type,
            model_name = 'model_true')

          # prepare starting model
          self.generate_mesh(
            model_path = PATH.MODEL_INIT,
            model_type = model_type,
            model_name = 'model_init')



  ### high-level solver interface

  def evaluate_func(self,path='',export_traces=False):
      unix.cd(self.getpath())
      self.import_model(path)

      # forward simulation
      self.forward()
      unix.mv(self.glob(),'traces/syn')
      self.preprocess.prepare_adjoint(unix.pwd(),output_type=2)

      # save results
      self.export_residuals(path)
      if export_traces:
        self.export_traces(path,'traces/syn')


  def evaluate_grad(self,path='',export_traces=False):
      unix.cd(self.getpath())

      # adjoint simulation
      self.adjoint()

      # save results
      self.export_kernels(path)
      if export_traces:
        self.export_traces(path,'traces/syn')


  def apply_hess(self,path='',hessian='exact'):
      unix.cd(self.getpath())
      self.import_model(path)

      # forward simulation
      self.forward()
      unix.mv(self.glob(),'traces/lcg')
      if PAR.SCHEME in ['gn']:
        self.preprocess.prepare_adjoint(unix.pwd(),output_type=3)
      elif PAR.SCHEME in ['tn']:
        self.preprocess.prepare_adjoint(unix.pwd(),output_type=4)

      # adjoint simulation
      self.adjoint()

      # save results
      self.export_kernels(path)



  ### low-level solver interface

  def forward(self):
      """ Calls forward solver
      """
      seistools.specfem2d.setpar('SIMULATION_TYPE', '1')
      seistools.specfem2d.setpar('SAVE_FORWARD',  ' .true.')

      f = open('log.fwd','w')
      subprocess.call(self.mesher_binary,stdout=f)
      subprocess.call(self.solver_binary,stdout=f)
      f.close()


  def adjoint(self):
      """ Calls adjoint solver
      """
      seistools.specfem2d.setpar('SIMULATION_TYPE', '3')
      seistools.specfem2d.setpar('SAVE_FORWARD',    '.false.')
 
      unix.rm('SEM')
      unix.ln('traces/adj','SEM')

      f = open('log.adj','w')
      subprocess.call(self.mesher_binary,stdout=f)
      subprocess.call(self.solver_binary,stdout=f)
      f.close()


  def generate_data(self,model_path='',model_type='',model_name=''):
      """ Generates data using SPECFEM2D forward solver
      """
      unix.cd(self.getpath())

      # prepare model
      self.generate_mesh(model_path,model_type,model_name)

      # prepare solver
      s = np.loadtxt('DATA/sources.dat')[system.getnode(),:]
      seistools.specfem2d.setpar('xs',s[0],file='DATA/SOURCE')
      seistools.specfem2d.setpar('zs',s[1],file='DATA/SOURCE')
      seistools.specfem2d.setpar('NSOURCES',str(1))

      seistools.specfem2d.setpar('SIMULATION_TYPE', '1')
      seistools.specfem2d.setpar('SAVE_FORWARD', ' .false.')

      # run solver
      f = open('log.fwd','w')
      subprocess.call(self.mesher_binary,stdout=f)
      subprocess.call(self.solver_binary,stdout=f)
      f.close()

      # export traces
      unix.mv(self.glob(),'traces/obs')
      self.exprt(PATH.OUTPUT,'obs')


  def generate_mesh(self,model_path='',model_type='',model_name=''):
      """ Performs meshing using SPECFEM2D builtin mesher and interpolates
          model parameters at mesh coordinates
      """
      assert exists(model_path)

      from copy import copy
      from time import sleep
      from scipy.interpolate import griddata
      import warnings
      warnings.filterwarnings('ignore')
      from seisflows.seistools.models import loadascii, loadnc, loadsep

      unix.cd(self.getpath())

      # load model
      if model_type in ['gll']:
        pass
      elif model_type in ['ascii','txt']:
        model = loadascii(model_path)
      elif model_type in ['h','h@']:
        model = loadsep(model_path)
      elif model_type in ['nc','netcdf','cdf','grd']:
        model = loadnc(model_path)
      else:
        raise Exception

      if model_type in ['gll']:
          parts = self.load(model_path)
          unix.cp(model_path,'DATA/model_velocity.dat_input')
      else:
        seistools.specfem2d.setpar('SIMULATION_TYPE', '1')
        seistools.specfem2d.setpar('SAVE_FORWARD', '.false.')
        seistools.specfem2d.setpar('assign_external_model', '.false.')
        seistools.specfem2d.setpar('READ_EXTERNAL_SEP_FILE', '.false.')
        nt = seistools.specfem2d.getpar('nt')
        seistools.specfem2d.setpar('nt',str(100))

        mesher = 'bin/xmeshfem2D'
        solver = 'bin/xspecfem2D'
        with open('log.adj','w') as f:
          subprocess.call(mesher,stdout=f)
          subprocess.call(solver,stdout=f)

        seistools.specfem2d.setpar('assign_external_model', '.true.')
        seistools.specfem2d.setpar('READ_EXTERNAL_SEP_FILE', '.true.')
        seistools.specfem2d.setpar('nt',str(nt))

        # interpolate model at GLL points
        parts = self.load('DATA/model_velocity.dat_output')
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

      # save results
      if system.getnode()==0:
        if model_name:
          self.save(PATH.OUTPUT+'/'+model_name,parts)
        if not exists(PATH.MESH):
          set1 = set(['x','z','rho','vp','vs'] )
          set2 = set(self.inversion_parameters)
          keys = list(set1.difference(set2))
          for key in keys:
            unix.mkdir(PATH.MESH+'/'+key)
            for proc in range(PAR.NPROC):
              with open(PATH.MESH+'/'+key+'/'+'%06d'%proc,'w') as file:
                np.save(file,parts[key][proc])



  ### model input/output

  def load(self,filename,type=''):
      "reads SPECFEM2D kernel or model"
      # read text file
      M = np.loadtxt(filename)
      nrow = M.shape[0]
      ncol = M.shape[1]

      if ncol==5:
        ioff = 0
      elif ncol==6:
        ioff = 1

      # fill in dictionary
      parts = {}
      for key in ['x','z','rho','vp','vs']:
        parts[key] = [M[:,ioff]]
        ioff += 1
      return parts


  def save(self,filename,parts,type='model'):
      "writes SPECFEM2D kernel or model"
      # allocate array
      if type == 'model':
        nrow = len(parts[parts.keys().pop()][0])
        ncol = 6
        ioff = 1
        M = np.zeros((nrow,ncol))

      elif type =='kernel':
        nrow = len(parts[parts.keys().pop()][0])
        ncol = 5
        ioff = 0
        M = np.zeros((nrow,ncol))

      # fill in array
      for icol,key in enumerate(['x','z','rho','vp','vs']):
        M[:,icol+ioff] = parts[key][0]

      # write array
      np.savetxt(filename,M,'%10.4e')



  ### vector/dictionary conversion

  def merge(self,parts):
      "merges parts into vector"
      v = np.array([])
      for key in parts:
        if key in self.inversion_parameters:
          for iproc in range(PAR.NPROC):
            v = np.append(v,parts[key][iproc])
      return v


  def split(self,v):
      "split vector into parts"
      parts = {}
      nrow = len(v)/(PAR.NPROC*len(self.inversion_parameters))
      i = 0; j = 0;
      for key in ['x','z','rho','vp','vs']:
        parts[key] = []
        if key in self.inversion_parameters:
          for i in range(PAR.NPROC):
            imin = nrow*PAR.NPROC*j + nrow*i
            imax = nrow*PAR.NPROC*j + nrow*(i+1)
            i += 1
            parts[key].append(v[imin:imax])
          j += 1
        else:
          for i in range(PAR.NPROC):
            proc = '%06d'%i
            parts[key].append(np.load(PATH.MESH+'/'+key+'/'+proc))
      return parts



  ### model manipulation

  def combine(self,path=''):
      "combines SPECFEM2D kernels"
      subprocess.call(\
             [getpath('seistools')+'/'+'specfem2d/combine.exe'] + \
             [str(len(unix.ls(path)))] + \
             [path])


  def smooth(self,path='',span=0):
      "smooths SPECFEM2D kernels"
      from seisflows.tools.arraytools import meshsmooth

      parts = self.load(path+'/'+'grad')
      if not span:
        return parts

      # set up grid
      x = parts['x'][0]
      z = parts['z'][0]
      lx = x.max() - x.min()
      lz = z.max() - z.min()
      nn = x.size
      nx = np.around(np.sqrt(nn*lx/lz))
      nz = np.around(np.sqrt(nn*lx/lz))

      # perform smoothing
      for key in self.inversion_parameters:
        parts[key] = [meshsmooth(x,z,parts[key][0],span,nx,nz)]
      unix.mv(path+'/'+'grad',path+'/'+'grad_nosmooth')
      self.save(path+'/'+'grad',parts)



  ### input file writers

  def writepar(self):
      "Write solver parameters file"
      seistools.specfem2d.writepar(PAR.vars)


  def writerec(self,prefix='traces/obs'):
      "Write receivers file"
      # adjust parameters
      key = 'use_existing_STATIONS'
      val = '.true.'
      seistools.specfem2d.setpar(key,val)

      # write receiver file
      _,h = self.preprocess.load(prefix)
      seistools.specfem2d.writerec(h.nr,h.rx,h.rz)


  def writesrc(self):
      "Write sources file"
      _,h = self.preprocess.load(prefix='traces/obs')
      seistools.specfem2d.writesrc(PAR.vars,h)



  ### file transfer utilities

  def import_model(self,path):
      src=join(path,'model')
      dst=join(unix.pwd(),'DATA/model_velocity.dat_input')
      unix.cp(src,dst)

  def import_traces(self,path):
      src=glob(join(path,'traces',getname(),'*'))
      dst=join(unix.pwd(),'traces/obs')
      unix.cp(src,dst)

  def export_kernels(self,path):
      try: unix.mkdir(join(path,'kernels'))
      except: pass
      src=join(unix.pwd(),'OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat')
      dst=join(path,'kernels','%06d'%system.getnode())
      unix.cp(src,dst)

  def export_residuals(self,path):
      try: unix.mkdir(join(path,'residuals'))
      except: pass
      src=join(unix.pwd(),'residuals')
      dst=join(path,'residuals','%06d'%system.getnode())
      unix.cp(src,dst)

  def export_traces(self,path,prefix='traces/obs'):
      try: unix.mkdir(join(path,'traces'))
      except: pass
      src=join(unix.pwd(),prefix)
      dst=join(path,'traces','%06d'%system.getnode())
      unix.cp(src,dst)



  ### utility functions

  def cleanup(self):
      """ Cleans up local directory after simulation
      """
      unix.cd(self.getpath())
      unix.rm(glob('traces/syn/*'))
      unix.rm(glob('traces/adj/*'))


  def getpath(self):
     return join(PATH.SCRATCH,self.getshot())


  def getshot(self):
    return '%06d'%system.getnode()



  ### class variables

  mesher_binary = 'bin/xmeshfem2D'
  solver_binary = 'bin/xspecfem2D'
  glob = lambda _ : glob('OUTPUT_FILES/U?_file_single.bin')

