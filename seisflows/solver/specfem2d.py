
import subprocess

import numpy as np
from scipy.interpolate import griddata

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join, setdiff, Struct
from seisflows.tools.configtools import loadclass, findpath, ParameterObj, ConfigObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class specfem2d(object):
    """ Python interface for SPECFEM2D

      evaluate_func, evaluate_grad, apply_hess
        These methods deal with evaluation of the misfit function or its
        derivatives and provide the primary interface between the solver and other
        workflow components.

      forward, adjoint, mesher
        These methods allow direct access to individual SPECFEM2D components.
        Together, they provide a secondary interface users can employ for
        specialized tasks not covered by high level methods.

      prepare_solver, prepare_data, prepare_model
        SPECFEM2D requires a particular directory structure in which to run and
        particular file formats for models, data, and parameter files. These
        methods help put in place all these prerequisites.

      load, save
        For reading and writing SPECFEM2D models and kernels. On the disk, models
        and kernels are stored as text files, and in memory, as dictionaries with
        different keys corresponding to different material parameters.

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
    inversion_parameters += ['vs']

    # data channels
    channels = []
    channels += ['y']

    # input/output
    reader = staticmethod(seistools.specfem2d.readsu)
    writer = staticmethod(seistools.specfem2d.writesu)

    mesher_binary = 'bin/xmeshfem2D'
    solver_binary = 'bin/xspecfem2D'

    glob = lambda _ : glob('OUTPUT_FILES/U?_file_single.bin')


    def check(self):

        global preprocess
        import preprocess

        global system
        import system

        # check scratch paths
        if 'GLOBAL' not in PATH:
            raise Exception

        if 'LOCAL' not in PATH:
            setattr(PATH,'LOCAL',None)

        if 'MESH' not in PATH:
            setattr(PATH,'MESH',join(PATH.GLOBAL,'mesh'))

        if 'SOLVER' not in PATH:
            if PATH.LOCAL:
                setattr(PATH,'SOLVER',join(PATH.LOCAL,'solver'))
            else:
                setattr(PATH,'SOLVER',join(PATH.GLOBAL,'solver'))

        # check mesh parameters
        if 'XMIN' not in PAR:
            raise Exception

        if 'XMAX' not in PAR:
            raise Exception

        if 'ZMIN' not in PAR:
            raise Exception

        if 'ZMAX' not in PAR:
            raise Exception

        if 'NX' not in PAR:
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


    def prepare_solver(self,inversion=True,**kwargs):
        """ Sets up directory from which to run solver

          Creates subdirectories expected by SPECFEM2D, copies mesher and solver
          binary files, and optionally calls prepare_model and prepare_data.
          Binaries must be supplied by user as there is currently no mechanism to
          automatically compile from source code.
        """
        unix.rm(self.path)
        unix.mkdir(self.path)
        unix.cd(self.path)

        # create subdirectories
        unix.mkdir('bin')
        unix.mkdir('DATA')
        unix.mkdir('OUTPUT_FILES')
        unix.mkdir('traces/obs')
        unix.mkdir('traces/syn')
        unix.mkdir('traces/adj')

        # copy binaries
        assert exists(PATH.SOLVER_BINARIES)
        src = glob(PATH.SOLVER_BINARIES+'/'+'*')
        dst = 'bin/'
        unix.cp(src,dst)

        if inversion:
            if 'model_type' not in kwargs:
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
          generated on the fly. If data are provided, SPECFEM2D input files will
          be generated automatically from included header information. If data
          are to be generated on the fly, both a target model and SPECFEM2D input
          files must be supplied.

          Adjoint traces are intialized by writing zeros for all components. This
          is necessary because, at the start of an adjoint simulation, SPECFEM2D
          expects all components to exist, even ones not actually in use for the
          inversion.
        """
        unix.cd(self.path)

        if exists(PATH.DATA):
            # copy user supplied data
            src = glob(PATH.DATA+'/'+self.getshot()+'/'+'*')
            dst = 'traces/obs/'
            unix.cp(src,dst)

            # generate SPECFEM2D input files
            self.write_parameters()
            self.write_sources()
            self.write_receivers()

        else:
            # generate data
            self.generate_data(**kwargs)

        # initialize adjoint traces
        zeros = np.zeros((PAR.NT,PAR.NREC))
        _,h = preprocess.load('traces/obs')
        for channel in ['x','y','z']:
            self.writer(zeros,h,channel=channel,prefix='traces/adj')


    def prepare_model(self,model_path='',model_type='',model_name=''):
        """ Performs meshing and model interpolation

          Meshing is carried out using SPECFEM2D's builtin mesher. Following that,
          model values are read from a user supplied input file. If the
          coordinates on which the model is defined do not match the coordinates
          of the mesh, an additional interpolation step is performed.
        """
        assert(exists(model_path))
        assert(model_type)

        unix.cd(self.path)

        # perform meshing
        if model_type in ['gll']:
            parts = self.load(model_path)
            unix.cp(model_path,'DATA/model_velocity.dat_input')
        else:
            self.mesher()

        # load material parameters
        if model_type in ['gll']:
            parts = self.load(model_path)
        else:
            if model_type in ['ascii','txt']:
                model = seistools.core.loadascii(model_path)
            elif model_type in ['h','h@']:
                model = seistools.core.loadsep(model_path)
            elif model_type in ['nc','netcdf','cdf','grd']:
                model = seistools.core.loadnc(model_path)
            else:
                raise Exception
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

        # save results
        if system.getnode()==0:
            if model_name:
                self.save(PATH.OUTPUT+'/'+model_name,parts)
            if not exists(PATH.MESH):
                for key in setdiff(['x','z','rho','vp','vs'],self.inversion_parameters):
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
        assert(exists(PATH.SOLVER_FILES))

        unix.cd(self.path)

        # copy input files
        src = glob(PATH.SOLVER_FILES+'/'+'*')
        dst = 'DATA/'
        unix.cp(src,dst)

        # set source coordinates
        s = np.loadtxt('DATA/sources.dat')[system.getnode(),:]
        seistools.specfem2d.setpar('xs',s[0],file='DATA/SOURCE')
        seistools.specfem2d.setpar('zs',s[1],file='DATA/SOURCE')
        seistools.specfem2d.setpar('NSOURCES',str(1))

        # prepare model
        self.prepare_model(**kwargs)

        # forward simulation
        seistools.specfem2d.setpar('SIMULATION_TYPE', '1')
        seistools.specfem2d.setpar('SAVE_FORWARD', '.false.')
        with open('log.fwd','w') as f:
            subprocess.call(self.mesher_binary,stdout=f)
            subprocess.call(self.solver_binary,stdout=f)

        # export traces
        unix.mv(self.glob(),'traces/obs')
        self.export_traces(PATH.OUTPUT,'traces/obs')



    ### high-level solver interface

    def evaluate_func(self,path='',export_traces=False):
        """ Evaluates misfit function by carrying out forward simulation and
            making measurements on observations and synthetics.
        """
        unix.cd(self.path)
        self.import_model(path)

        # forward simulation
        self.forward()
        unix.mv(self.glob(),'traces/syn')
        preprocess.prepare_adjoint(unix.pwd(),output_type=2)

        # save results
        self.export_residuals(path)
        if export_traces:
            self.export_traces(path,'traces/syn')


    def evaluate_grad(self,path='',export_traces=False):
        """ Evaluates gradient by carrying out adjoint simulation. Adjoint traces
            must already be in place prior to calling this method or the adjoint
            simulation will fail.
        """
        unix.cd(self.path)

        # adjoint simulation
        self.adjoint()

        # save results
        self.export_kernels(path)
        if export_traces:
            self.export_traces(path,'traces/syn')


    def apply_hess(self,path='',hessian='Newton'):
        """ Evaluates action of Hessian on a given model vector.
        """
        unix.cd(self.path)
        self.import_model(path)

        # forward simulation
        self.forward()
        unix.mkdir('traces/lcg')
        unix.mv(self.glob(),'traces/lcg')

        if hessian == 'Newton':
            preprocess.prepare_adjoint(unix.pwd(),output_type=3)

        elif hessian == 'GaussNewton':
            preprocess.prepare_adjoint(unix.pwd(),output_type=4)

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

        with open('log.fwd','w') as f:
            subprocess.call(self.mesher_binary,stdout=f)
            subprocess.call(self.solver_binary,stdout=f)


    def adjoint(self):
        """ Calls adjoint solver
        """
        seistools.specfem2d.setpar('SIMULATION_TYPE', '3')
        seistools.specfem2d.setpar('SAVE_FORWARD',    '.false.')

        unix.rm('SEM')
        unix.ln('traces/adj','SEM')

        with  open('log.adj','w') as f:
            subprocess.call(self.mesher_binary,stdout=f)
            subprocess.call(self.solver_binary,stdout=f)


    def mesher(self):
        """ Calls SPECFEM2D builtin mesher
        """
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


    ### model input/output

    def load(self,filename,type=''):
        """Reads SPECFEM2D kernel or model

           Models and kernels are read from 5 or 6 column text files whose format
           is described in the SPECFEM2D user manual. Once read, a model or kernel
           is stored in a dictionary containing mesh coordinates and
           corresponding material parameter values.
        """
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



    ### postprocessing utilities

    def combine(self,path=''):
        "combines SPECFEM2D kernels"
        subprocess.call(\
               [findpath('seistools')+'/'+'specfem2d/combine.exe'] + \
               [str(len(unix.ls(path)))] + \
               [path])


    def smooth(self,path='',tag='grad',span=0):
        "smooths SPECFEM2D kernels by convoling them with a Gaussian"
        from seisflows.tools.arraytools import meshsmooth

        parts = self.load(path+'/'+tag)
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
        unix.mv(path+'/'+tag,path+'/'+tag+'_nosmooth')
        self.save(path+'/'+tag,parts)



    ### input file writers

    def write_parameters(self):
        "Write solver parameters file"
        unix.cd(self.path)
        seistools.specfem2d.write_parameters(vars(PAR))


    def write_receivers(self,prefix='traces/obs'):
        "Write receivers file"
        unix.cd(self.path)

        # adjust parameters
        key = 'use_existing_STATIONS'
        val = '.true.'
        seistools.specfem2d.setpar(key,val)

        # write receiver file
        _,h = preprocess.load(prefix)
        seistools.specfem2d.write_receivers(h.nr,h.rx,h.rz)


    def write_sources(self):
        "Write sources file"
        unix.cd(self.path)

        _,h = preprocess.load(prefix='traces/obs')
        seistools.specfem2d.write_sources(vars(PAR),h)



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

    @property
    def path(self):
         return join(PATH.SOLVER,self.getshot())


    def getshot(self):
        return '%06d'%system.getnode()

