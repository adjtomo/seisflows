
import subprocess

import numpy as np

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configure import getclass, getpath, ParameterObject
from seisflows.tools.iotools import loadbin, savebin


PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class specfem3d(object):
    """ Python interface for SPECFEM3D

      evaluate_func, evaluate_grad, apply_hess
        These methods deal with evaluation of the misfit function or its
        derivatives and provide the primary interface between the solver and other
        workflow components.

      forward, adjoint, mesher
        These methods allow direct access to individual SPECFEM3D components.
        Together, they provide a secondary interface users can employ for
        specialized tasks not covered by high level methods.

      prepare_solver, prepare_data, prepare_model
        SPECFEM3D requires a particular directory structure in which to run and
        particular file formats for models, data, and parameter files. These
        methods help put in place all these prerequisites.

      load, save
        For reading and writing SPECFEM3D models and kernels. On the disk, models
        and kernels are stored as binary files, and in memory, as dictionaries
        with different keys corresponding to different material parameters.

      split, merge
        In the solver routines, it is possible to store models as dictionaries,
        but for the optimization routines, it is necessary to merge all model
        values together into a single vector. Two methods, 'split' and 'merge',
        are used to convert back and forth between these two representations.

      combine, smooth
        Utilities for combining and smoothing kernels, meant to be called from
        external postprocessing routines.
    """

    # model parameters
    model_parameters = []
    model_parameters += ['rho']
    model_parameters += ['vp']
    model_parameters += ['vs']

    # inversion parameters
    inversion_parameters = []
    inversion_parameters += ['vp']
    inversion_parameters += ['vs']

    # data channels
    channels = []
    channels += ['z']

    kernel_map = {
      'rho':'rho_kernel',
      'vp':'alpha_kernel',
      'vs':'beta_kernel'}

    def __init__(self):
        """ Checks user supplied parameters
        """
        # check mesh parameters
        if 'MESH' not in PATH:
            setattr(PATH,'MESH',join(PATH.SOLVER,'mesh'))

        if 'XMIN' not in PAR:
            raise Exception

        if 'XMAX' not in PAR:
            raise Exception

        if 'YMIN' not in PAR:
            raise Exception

        if 'YMAX' not in PAR:
            raise Exception

        if 'ZMIN' not in PAR:
            raise Exception

        if 'ZMAX' not in PAR:
            raise Exception

        if 'NX' not in PAR:
            raise Exception

        if 'NY' not in PAR:
            raise Exception

        if 'NZ' not in PAR:
            raise Exception

        # check time stepping parameters
        if 'NT' not in PAR:
            raise Exception

        if 'DT' not in PAR:
            raise Exception

        if 'F0' not in PAR:
            raise Exception

        if 'WAVELET' not in PAR:
            setattr(PAR,'WAVELET','ricker')

        # load preprocessing machinery
        if 'PREPROCESS' not in PAR:
            setattr(PAR,'PREPROCESS','default')

        self.preprocess = getclass('preprocess',PAR.PREPROCESS)(
          channels=self.channels,
          reader=seistools.specfem3d.readsu,
          writer=seistools.specfem3d.writesu)


    def prepare_solver(self,inversion=True,**kwargs):
        """ Sets up directory from which to run solver

          Creates subdirectories expected by SPECFEM3D, copies mesher and solver
          binary files, and optionally calls prepare_data and prepare_mdoel.
          Binaries must be supplied by user as there is currently no mechanism to
          automatically compile from source code.
        """
        unix.rm(self.path())
        unix.mkdir(self.path())
        unix.cd(self.path())

        # create subdirectories
        unix.mkdir('bin')
        unix.mkdir('DATA')
        unix.mkdir('OUTPUT_FILES/DATABASES_MPI')

        unix.mkdir('traces/obs')
        unix.mkdir('traces/syn')
        unix.mkdir('traces/adj')

        # copy binaries
        src = glob(PATH.SOLVER_BINARIES+'/'+'*')
        dst = 'bin/'
        unix.cp(src,dst)

        if inversion:
            if 'model_type' not in kwargs:
                try:
                    file = PATH.SOLVER_FILES+'/'+'Par_file'
                    kwargs['model_type'] = seistools.specfem3d.getpar('MODEL',file=file)
                except:
                    kwargs['model_type'] = 'gll'

            # prepare data
            if not exists(PATH.DATA):
                assert(exists(PATH.SOLVER_FILES))
                kwargs['model_path'] = PATH.MODEL_TRUE
                kwargs['model_name'] = 'model_true'
            self.prepare_data(**kwargs)

            # prepare model
            kwargs['model_path'] = PATH.MODEL_INIT
            kwargs['model_name'] = 'model_init'
            self.prepare_model(**kwargs)


    def prepare_data(self,**kwargs):
        """ Prepares data for inversion or migration

          Users implementing an inversion or migration can choose between
          supplying data or supplying a target model from which data are
          generated on the fly. If data are provided, SPECFEM3D input files will
          be generated automatically from included header information. If data
          are to be generated on the fly, both a target model and SPECFEM3D input
          files must be supplied.

          Adjoint traces are intialized by writing zeros for all components. This
          is necessary because, at the start of an adjoint simulation, SPECFEM3D
          expects that all components exists, even ones not actually in use for
          the inversion.
        """
        unix.cd(self.path())

        if exists(PATH.DATA):
                # copy user supplied data
            src = glob(PATH.DATA+'/'+getname()+'/'+'*')
            dst = 'traces/obs/'
            unix.cp(src,dst)

            # generate SPECFEM3D input files
            self.write_parameters()
            self.write_sources()
            self.write_receivers()

        else:
            # generate data
            self.generate_data(**kwargs)

        # initialize adjoint traces
        zeros = np.zeros((PAR.NT,PAR.NREC))
        _,h = self.preprocess.load('traces/obs')
        for channel in ['x','y','z']:
            self.preprocess.writer(zeros,h,
                channel=channel,prefix='traces/adj')


    def prepare_model(self,model_path='',model_type='',model_name=''):
        """ Performs meshing and model interpolation using SPECFEM3D's builtin
          mesher and database generation utility
        """
        assert(model_type)

        unix.cd(self.path())

        # run builtin mesher and generate databases
        if model_type == 'gll':
            assert(exists(model_path))
            # copy files
            src = glob(model_path+'/'+'*')
            dst = 'OUTPUT_FILES/DATABASES_MPI/'
            unix.cp(src,dst)

        elif model_type == 'sep':
            assert(exists(model_path))
            # copy files
            src = glob(model_path+'/'+'*')
            dst = 'DATA/'
            unix.cp(glob(model_path+'/'+'*'),'DATA/')

        elif model_type == 'default':
            pass

        elif model_type == 'tomo':
            pass

        seistools.specfem3d.setpar('MODEL',model_type)
        self.mesher()

        # save results
        parts = self.load('OUTPUT_FILES/DATABASES_MPI')
        if system.getnode()==0:
            if model_name and model_type == 'gll':
                unix.ln(model_path,PATH.OUTPUT+'/'+model_name)
            elif model_name:
                self.save(PATH.OUTPUT+'/'+model_name,parts)
            if not exists(PATH.MESH):
                set1 = set(self.model_parameters)
                set2 = set(self.inversion_parameters)
                keys = list(set1.difference(set2))
                for key in keys:
                    unix.mkdir(PATH.MESH+'/'+key)
                    for proc in range(PAR.NPROC):
                        with open(PATH.MESH+'/'+key+'/'+'%06d'%proc,'w') as file:
                            np.save(file,parts[key][proc])


    def generate_data(self,**kwargs):
        """ Generates data by copying input files, setting source coordinates,
          preparing model, and  calling forward solver. Similar to 'forward',
          but puts together all the steps needed for modeling, rather than just
          calling the solver binary.
        """
        unix.cd(self.path())

        # copy input files
        src = glob(PATH.SOLVER_FILES+'/'+'*')
        dst = 'DATA/'
        unix.cp(src,dst)

        # set source coordinates
        s = np.loadtxt('DATA/SOURCE.XYZ')[system.getnode(),:]
        seistools.specfem3d.setpar('longitude',s[0],sep=': ',file='DATA/FORCESOLUTION')
        seistools.specfem3d.setpar('latitude',s[1],sep=': ',file='DATA/FORCESOLUTION')
        seistools.specfem3d.setpar('depth',s[2],sep=': ',file='DATA/FORCESOLUTION')

        # prepare model
        self.prepare_model(**kwargs)

        # run solver
        seistools.specfem3d.setpar('SIMULATION_TYPE', '1')
        seistools.specfem3d.setpar('SAVE_FORWARD', '.false.')
        self.mpirun('bin/xspecfem3D')
        unix.mv(self.glob(),'traces/obs')



    ### high-level solver interface

    def evaluate_func(self,path='',export_traces=False):
        """ Evaluates misfit function by carrying out forward simulation and
            making measurements on observations and synthetics.
        """
        unix.cd(self.path())
        self.import_model(path)

        # forward simulation
        self.forward()
        unix.mv(self.glob(),'traces/syn')
        self.preprocess.prepare_adjoint(unix.pwd(),output_type=2)

        # save results
        self.export_residuals(path)
        if export_traces:
            self.export_traces(path,prefix='traces/syn')


    def evaluate_grad(self,path='',export_traces=False):
        """ Evaluates gradient by carrying out adjoint simulation. Adjoint traces
            must already be in place prior to calling this method or the adjoint
            simulation will fail.
        """
        unix.cd(self.path())

        # adjoint simulation
        self.adjoint()

        # save results
        self.export_kernels(path)
        if export_traces:
            self.export_traces(path,prefix='traces/syn')


    def apply_hess(self,path='',hessian='exact'):
        """ Evaluates action of Hessian on a given model vector.
        """
        unix.cd(self.path())
        self.imprt(path,'model')

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
        self.export_kernels(path,'kernels')



    ### low-level solver interface

    def forward(self):
        """ Calls SPECFEM3D forward solver
        """
        # prepare solver
        seistools.specfem3d.setpar('SIMULATION_TYPE', '1')
        seistools.specfem3d.setpar('SAVE_FORWARD', '.true.')
        seistools.specfem3d.setpar('MODEL','gll')
        self.mpirun('bin/xgenerate_databases')

        # run solver
        self.mpirun('bin/xspecfem3D')
        unix.mv(self.glob(),'traces/syn')


    def adjoint(self):
        """ Calls SPECFEM3D adjoint solver
        """
        # prepare solver
        seistools.specfem3d.setpar('SIMULATION_TYPE', '3')
        seistools.specfem3d.setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj','SEM')

        # run solver
        self.mpirun('bin/xspecfem3D')


    def mesher(self):
        """ Calls SPECFEM3D builtin mesher
        """
        self.mpirun('bin/xmeshfem3D')
        self.mpirun('bin/xgenerate_databases')



    ### model input/output

    def load(self,dirname,type='model'):
        """ reads SPECFEM3D kernel or model to dictionary
        """
        mapping = lambda key : self.kernel_map[key]
        parts = {}

        if type == 'model':
            for key in self.model_parameters:
                parts[key] = []
                # read database files
                for iproc in range(PAR.NPROC):
                    filename = 'proc%06d_%s.bin' % (iproc,key)
                    part = loadbin(join(dirname,filename))
                    parts[key].append(part)

        elif type == 'kernel':
            for key in self.model_parameters:
                parts[key] = []
                # read database files
                for iproc in range(PAR.NPROC):
                    filename = 'proc%06d_%s.bin' % (iproc,mapping(key))
                    part = loadbin(join(dirname,filename))
                    parts[key].append(part)

        return parts


    def save(self,dirname,parts):
        """ writes SPECFEM3D model
        """
        unix.mkdir(dirname)

        # write database files
        for key in self.model_parameters:
            nn = len(parts[key])
            for ii in range(nn):
                filename = 'proc%06d_%s.bin' % (ii,key)
                savebin(parts[key][ii],join(dirname,filename))



    ### vector/dictionary conversion

    def merge(self,parts):
        """ merges dictionary into vector
        """
        v = np.array([])
        for key in self.inversion_parameters:
            for iproc in range(PAR.NPROC):
                v = np.append(v,parts[key][iproc])
        return v


    def split(self,v):
        """ splits vector into dictionary
        """
        parts = {}
        nrow = len(v)/(PAR.NPROC*len(self.inversion_parameters))
        i = 0; j = 0
        for key in self.model_parameters:
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



    ### postprocessing utilities

    def combine(self,path=''):
        """ combines SPECFEM3D kernels
        """
        unix.cd(self.path())

        # create temporary files and directories needed by xsum_kernels
        dirs = unix.ls(path)
        with open('../kernels_list.txt','w') as file:
            file.write('\n'.join(dirs)+'\n')
        unix.mkdir('INPUT_KERNELS')
        unix.mkdir('OUTPUT_SUM')
        for dir in dirs:
            src = path + '/' + dir
            dst = unix.pwd() +'/'+ 'INPUT_KERNELS' +'/'+ dir
            unix.ln(src,dst)

        # sum kernels
        self.mpirun(PATH.SOLVER_BINARIES+'/'+'xsum_kernels')
        unix.mv('OUTPUT_SUM',path+'/'+'sum')

        # remove temporary files and directories
        unix.rm('INPUT_KERNELS')
        unix.rm('../kernels_list.txt')

        unix.cd(path)


    def smooth(self,path='',span=0):
        """ smooths SPECFEM3D kernels
        """
        unix.cd(self.path())

        unix.mv(path+'/'+'gradient',path+'/'+'gradient_nosmooth')
        unix.mkdir(path+'/'+'gradient')

        # prepare list
        kernel_list = []
        for key in self.model_parameters:
            if key in self.inversion_parameters:
                smoothing_flag = True
            else:
                smoothing_flag = False
            kernel_list = kernel_list + [[key,smoothing_flag]]

        for kernel_name,smoothing_flag in kernel_list:
            if smoothing_flag:
                # run smoothing
                print ' smoothing', kernel_name
                self.mpirun(
                  PATH.SOLVER_BINARIES+'/'+'xsmooth_vol_data '
                  + kernel_name + ' '
                  + path + '/' + 'gradient_nosmooth/' + ' '
                  + path + '/' + 'gradient/' + ' '
                  + str(span) + ' '
                  + str(span) + ' ' + '1')
            else:
                src = glob(path+'/'+'gradient_nosmooth/*'+kernel_name+'.bin')
                dst = path+'/'+'gradient/'
                unix.cp(src,dst)

        unix.rename('_smooth','',glob(path+'/'+'gradient/*_smooth.bin'))
        print ''

        unix.cd(path)



    ### input file writers

    def write_parameters(self):
        unix.cd(self.path())

        write_parameters(vars(PAR))


    def write_receivers(self):
        unix.cd(self.path())

        # adjust parameters
        key = 'use_existing_STATIONS'
        val = '.true.'
        seistools.specfem3d.setpar(key,val)

        # write receivers file
        _,h = self.preprocess.load('traces/obs')
        write_receivers(h.nr,h.rx,h.rz)


    def write_sources(self):
        unix.cd(self.path())

        # write source file
        _,h = self.preprocess.load(dir='traces/obs')
        write_sources(vars(PAR),h)



    ### file transfer utilities

    def import_model(self,path):
        src=glob(join(path,'model','*'))
        dst=join(unix.pwd(),'OUTPUT_FILES/DATABASES_MPI')
        unix.cp(src,dst)

    def import_traces(self,path):
        src=glob(join(path,'traces',getname(),'*'))
        dst=join(unix.pwd(),'traces/obs')
        unix.cp(src,dst)

    def export_model(self,path):
        if system.getnode()!=0:
            return
        for name in self.model_parameters:
            src=glob(join(unix.pwd(),'OUTPUT_FILES/DATABASES_MPI','*_'+name+'.bin'))
            dst=path
            unix.mkdir(dst)
            unix.cp(src,dst)

    def export_kernels(self,path):
        try: unix.mkdir(join(path,'kernels'))
        except: pass
        unix.mkdir(join(path,'kernels','%06d'%system.getnode()))
        for name in self.kernel_map.values():
            src=join(glob(unix.pwd()+'/'+'OUTPUT_FILES/DATABASES_MPI'+'/'+'*'+name+'.bin'))
            dst=join(path,'kernels','%06d'%system.getnode())
            unix.mv(src,dst)
        try:
            name = 'rhop_kernel'
            src=join(glob(unix.pwd()+'/'+'OUTPUT_FILES/DATABASES_MPI'+'/'+'*'+name+'.bin'))
            dst=join(path,'kernels','%06d'%system.getnode())
            unix.mv(src,dst)
        except:
            pass

    def export_residuals(self,path):
        try: unix.mkdir(join(path,'residuals'))
        except: pass
        src=join(unix.pwd(),'residuals')
        dst=join(path,'residuals','%06d'%system.getnode())
        unix.mv(src,dst)

    def export_traces(self,path,prefix='traces/obs'):
        try: unix.mkdir(join(path,'traces'))
        except: pass
        src=join(unix.pwd(),prefix)
        dst=join(path,'traces','%06d'%system.getnode())
        unix.cp(src,dst)



    ### utility functions

    def cleanup(self):
        """ Cleans up directory after simulation
        """
        unix.cd(self.path())
        unix.rm(glob('traces/syn/*'))
        unix.rm(glob('traces/adj/*'))


    def mpirun(self,script,outfile='/dev/null'):
        """ Wrapper for mpirun
        """
        with open(outfile) as f:
            subprocess.call(
                  system.mpiargs()
                  + script,
                  shell=True,
                  stdout=f)


    def path(self):
        return join(PATH.SCRATCH,self.getshot())


    def getshot(self):
        return '%06d'%system.getnode()


    def gettype(self):
        try:
            return seistools.specfem3d.getpar(
              'MODEL',file=PATH.SOLVER_FILES+'/'+'Par_file')
        except:
            return None


    glob = lambda _ : glob('OUTPUT_FILES/*_SU')
