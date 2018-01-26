
import subprocess
import sys
import numpy as np

from functools import partial
from glob import glob
from importlib import import_module
from os.path import basename, join
from seisflows.config import ParameterError, custom_import
from seisflows.plugins import solver_io
from seisflows.tools import msg, unix
from seisflows.tools.seismic import Container, call_solver
from seisflows.tools.tools import Struct, diff, exists



PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']



class base(object):
    """ Provides an interface through which solver simulations can be set up
      and run and a parent class for SPECFEM2D, SPECFEM3D and SPECFEM3D_GLOBE 
      subclasses

      This class supports only acoustic and isotropic elastic inversions.
      For additional options, see github.com/rmodrak/seisflows-multiparameter

      eval_func, eval_grad, apply_hess
        These methods deal with evaluation of the misfit function or its
        derivatives.  Together, they provide the primary interface through which
        SeisFlows interacts with SPECFEM2D/3D

      forward, adjoint
        These methods allow direct access to low-level SPECFEM2D/3D components,
        providing an alternative interface through which to interact with the 
        solver

     steup, generate_data, generate_model
        One-time operations performed at beginning of an inversion or 
        migration

     initialize_solver_directories, initialize_adjoint_traces
        SPECFEM2D/3D requires a particular directory structure in which to run and
        particular file formats for models, data, and parameter files. These
        methods help put in place all these prerequisites

      load, save
        For reading and writing SPECFEM2D/3D models and kernels. On the disk,
        models and kernels are stored as binary files, and in memory, as
        dictionaries with different keys corresponding to different material
        parameters

      split, merge
        Within the solver routines, it is natural to store models as 
        dictionaries. Within the optimization routines, it is natural to store
        models as vectors. Two methods, 'split' and 'merge', are used to convert 
        back and forth between these two representations

      combine, smooth
        Utilities for combining and smoothing kernels

    """

    assert 'MATERIALS' in PAR
    assert 'DENSITY' in PAR

    parameters = []
    if PAR.MATERIALS == 'Elastic':
        parameters += ['vp']
        parameters += ['vs']
    elif PAR.MATERIALS == 'Acoustic':
        parameters += ['vp']

    if PAR.DENSITY == 'Variable':
        parameters += ['rho']


    def check(self):
        """ Checks parameters and paths
        """
        # number of processors per simulation
        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')


        # format used by SPECFEM for reading and writing models
        # (currently, SPECFEM offers both 'fortran_binary' and 'adios')
        if 'SOLVERIO' not in PAR:
            setattr(PAR, 'SOLVERIO', 'fortran_binary')


        # solver scratch paths
        if 'SCRATCH' not in PATH:
            raise ParameterError(PATH, 'SCRATCH')

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SOLVER' not in PATH:
            if PATH.LOCAL:
                setattr(PATH, 'SOLVER', join(PATH.LOCAL, 'solver'))
            else:
                setattr(PATH, 'SOLVER', join(PATH.SCRATCH, 'solver'))

        # solver input paths
        if 'SPECFEM_BIN' not in PATH:
            raise ParameterError(PATH, 'SPECFEM_BIN')

        if 'SPECFEM_DATA' not in PATH:
            raise ParameterError(PATH, 'SPECFEM_DATA')

        # assertions
        assert self.parameters != []
        assert hasattr(solver_io, PAR.SOLVERIO)
        assert hasattr(self.io, 'read_slice')
        assert hasattr(self.io, 'write_slice')


    def setup(self):
        """ 
          Prepares solver for inversion or migration
          Sets up directory structure expected by SPECFEM and copies or 
          generates seismic data to be inverted or migrated
        """
        # clean up for new inversion
        unix.rm(self.cwd)

        # As input for an inversion or migration, users can choose between
        # providing data, or providing a target model from which data are
        # generated on the fly. In the former case, a value for PATH.DATA must
        # be supplied; in the latter case, a value for PATH.MODEL_TRUE must be
        # provided

        if PATH.DATA:
            # copy user supplied data
            self.initialize_solver_directories()

            src = glob(PATH.DATA +'/'+ self.source_name +'/'+ '*')
            dst = 'traces/obs/'
            unix.cp(src, dst)

        else:
            # generate data on the fly
            self.generate_data(
                model_path=PATH.MODEL_TRUE,
                model_name='model_true',
                model_type='gll')

        # prepare initial model
        self.generate_mesh(
            model_path=PATH.MODEL_INIT,
            model_name='model_init',
            model_type='gll')

        self.initialize_adjoint_traces()


    def clean(self):
        unix.cd(self.cwd)
        unix.rm('OUTPUT_FILES')
        unix.mkdir('OUTPUT_FILES')


    def generate_data(self, *args, **kwargs):
        """ Generates data
        """
        # must be implemented by subclass
        raise NotImplementedError


    def generate_mesh(self, *args, **kwargs):
        """ Performs meshing and database generation
        """
        # must be implemented by subclass
        raise NotImplementedError


    ### high-level solver interface

    def eval_func(self, path='', export_traces=False, write_residuals=True):
        """
          Performs forward simulations needed for misfit function evaluation

          :input path :: directory from which model is imported
          :input export_traces :: save or discard traces?
        """
        unix.cd(self.cwd)
        self.import_model(path)
        self.forward()

        if write_residuals:
            preprocess.prepare_eval_grad(self.cwd)
            self.export_residuals(path)


    def eval_grad(self, path='', export_traces=False):
        """ 
          Evaluates gradient by carrying out adjoint simulations.
          (A function evaluation must already have been carried out.)

          :input path :: directory from which model is imported
          :input export_traces :: save or discard traces?
        """
        unix.cd(self.cwd)
        self.adjoint()
        self.export_kernels(path)
        if export_traces:
            self.export_traces(path+'/'+'traces/syn', prefix='traces/syn')
            self.export_traces(path+'/'+'traces/adj', prefix='traces/adj')


    def apply_hess(self, path=''):
        """
          Computes action of Hessian on a given model vector.
          (A gradient evaluation must have already been carried out.)
 
          :input path :: directory to which output files are exported
        """
        unix.cd(self.cwd)
        self.import_model(path)
        unix.mkdir('traces/lcg')
        self.forward('traces/lcg')
        preprocess.prepare_apply_hess(self.cwd)
        self.adjoint()
        self.export_kernels(path)



    ### low-level solver interface

    def forward(self):
        """ Calls forward solver
        """
        # must be implemented by subclass
        raise NotImplementedError


    def adjoint(self):
        """ Calls adjoint solver
        """
        # must be implemented by subclass
        raise NotImplementedError


    ### model input/output

    @property
    def io(self):
        """ Solver IO module
        """
        return getattr(solver_io, PAR.SOLVERIO)


    def load(self, path, parameters=[], prefix='', suffix=''):
        """ 
          Loads SPECFEM2D/3D models or kernels

          :input path :: directory from which model is read
          :input parameters :: list of material parameters to be read
              (if empty, defaults to self.parameters)
          :input prefix :: optional filename prefix
          :input suffix :: optional filename suffix, eg '_kernel'
          :output dict :: model or kernels indexed by material parameter
              and processor rank, ie dict[parameter][iproc]
        """
        dict = Container()
        for iproc in range(self.mesh_properties.nproc):
            for key in parameters or self.parameters:
                dict[key] += self.io.read_slice(
                    path, prefix+key+suffix, iproc)
        return dict


    def save(self, dict, path, parameters=['vp','vs','rho'],
             prefix='', suffix=''):
        """ 
          Saves SPECFEM2D/3D models or kernels

          :input dict :: model stored as a dictionary or Container
          :input path :: directory to which model is written
          :input parameters :: list of material parameters to be written
          :input prefix :: optional filename prefix
          :input suffix :: optional filename suffix, eg '_kernel'
        """
        unix.mkdir(path)

        # fill in any missing parameters
        missing_keys = diff(parameters, dict.keys())
        for iproc in range(self.mesh_properties.nproc):
            for key in missing_keys:
                dict[key] += self.io.read_slice(
                    PATH.MODEL_INIT, prefix+key+suffix, iproc)

        # write slices to disk
        for iproc in range(self.mesh_properties.nproc):
            for key in parameters:
                self.io.write_slice(
                    dict[key][iproc], path, prefix+key+suffix, iproc)


    def merge(self, model, parameters=[]):
        """ Converts model from dictionary to vector representation
        """
        m = np.array([])
        for key in parameters or self.parameters:
            for iproc in range(self.mesh_properties.nproc):
                m = np.append(m, model[key][iproc])
        return m


    def split(self, m, parameters=[]):
        """ Converts model from vector to dictionary representation
        """
        nproc = self.mesh_properties.nproc
        ngll = self.mesh_properties.ngll
        model = Container()
        for idim, key in enumerate(parameters or self.parameters):
            model[key] = []
            for iproc in range(nproc):
                imin = sum(ngll)*idim + sum(ngll[:iproc])
                imax = sum(ngll)*idim + sum(ngll[:iproc+1])
                model[key] += [m[imin:imax]]
        return model



    ### postprocessing wrappers

    def combine(self, input_path='', output_path='', parameters=[]):
        """ Sums individual source contributions. Wrapper over xcombine_sem
            utility.
        """
        if not exists(input_path):
            raise Exception

        if not exists(output_path):
            unix.mkdir(output_path)

        unix.cd(self.cwd)
        with open('kernel_paths', 'w') as file:
            file.writelines([join(input_path, name+'\n')
                for name in self.source_names])

        for name in parameters or self.parameters:
            call_solver(
                system.mpiexec(),
                PATH.SPECFEM_BIN +'/'+ 'xcombine_sem '
                + name + '_kernel' + ' '
                + 'kernel_paths' + ' '
                + output_path)


    def smooth(self, input_path='', output_path='', parameters=[], span=0.):
        """ Smooths kernels by convolving them with a Gaussian.  Wrapper over 
            xsmooth_sem utility.
        """
        if not exists(input_path):
            raise Exception

        if not exists(output_path):
            unix.mkdir(output_path)

        # apply smoothing operator
        unix.cd(self.cwd)
        for name in parameters or self.parameters:
            print ' smoothing', name
            call_solver(
                system.mpiexec(),
                PATH.SPECFEM_BIN +'/'+ 'xsmooth_sem '
                + str(span) + ' '
                + str(span) + ' '
                + name + '_kernel' + ' '
                + input_path + '/ '
                + output_path + '/ ',
                output='/dev/null')

        print ''

        # rename output files
        files = glob(output_path+'/*')
        unix.rename('_smooth', '', files)


    ### file transfer utilities

    def import_model(self, path):
        model = self.load(path+'/'+'model')
        self.save(model, self.model_databases)

    def import_traces(self, path):
        src = glob(join(path, 'traces', self.source_name, '*'))
        dst = join(self.cwd, 'traces/obs')
        unix.cp(src, dst)

    def export_model(self, path, parameters=['rho', 'vp', 'vs']):
        if self.taskid == 0:
            unix.mkdir(path)
            for key in parameters:
                files = glob(join(self.model_databases, '*'+key+'.bin'))
                unix.cp(files, path)

    def export_kernels(self, path):
        unix.cd(self.kernel_databases)

        # work around conflicting name conventions
        self.rename_kernels()

        src = glob('*_kernel.bin')
        dst = join(path, 'kernels', self.source_name)
        unix.mkdir(dst)
        unix.mv(src, dst)

    def export_residuals(self, path):
        unix.mkdir(join(path, 'residuals'))

        src = join(self.cwd, 'residuals')
        dst = join(path, 'residuals', self.source_name)
        unix.mv(src, dst)

    def export_traces(self, path, prefix='traces/obs'):
        unix.mkdir(join(path))

        src = join(self.cwd, prefix)
        dst = join(path, self.source_name)
        unix.cp(src, dst)


    def rename_kernels(self):
        """ Works around conflicting kernel filename conventions
        """
        files = []
        files += glob('*proc??????_alpha_kernel.bin')
        files += glob('*proc??????_alpha[hv]_kernel.bin')
        files += glob('*proc??????_reg1_alpha_kernel.bin')
        files += glob('*proc??????_reg1_alpha[hv]_kernel.bin')
        unix.rename('alpha', 'vp', files)

        files = []
        files += glob('*proc??????_beta_kernel.bin')
        files += glob('*proc??????_beta[hv]_kernel.bin')
        files += glob('*proc??????_reg1_beta_kernel.bin')
        files += glob('*proc??????_reg1_beta[hv]_kernel.bin')
        unix.rename('beta', 'vs', files)


    def rename_data(self, path):
        """ Works around conflicting data filename conventions
        """
        pass


    ### setup utilities

    def initialize_solver_directories(self):
        """ Creates directory structure expected by SPECFEM3D, copies 
          executables, and prepares input files. Executables must be supplied 
          by user as there is currently no mechanism for automatically
          compiling from source.
        """
        unix.mkdir(self.cwd)
        unix.cd(self.cwd)

        # create directory structure
        unix.mkdir('bin')
        unix.mkdir('DATA')
        unix.mkdir('OUTPUT_FILES')

        unix.mkdir('traces/obs')
        unix.mkdir('traces/syn')
        unix.mkdir('traces/adj')

        unix.mkdir(self.model_databases)
        unix.mkdir(self.kernel_databases)

        # copy exectuables
        src = glob(PATH.SPECFEM_BIN +'/'+ '*')
        dst = 'bin/'
        unix.cp(src, dst)

        # copy input files
        src = glob(PATH.SPECFEM_DATA +'/'+ '*')
        dst = 'DATA/'
        unix.cp(src, dst)

        src = 'DATA/' + self.source_prefix +'_'+ self.source_name
        dst = 'DATA/' + self.source_prefix
        unix.cp(src, dst)

        self.check_solver_parameter_files()


    def initialize_adjoint_traces(self):
        """ Puts in place "adjoint traces" expected by SPECFEM
        """
        for filename in self.data_filenames:
            # read traces
            d = preprocess.reader(self.cwd +'/'+ 'traces/obs', filename)

            # Adjoint traces are initialized by writing zeros for all channels.
            # Channels actually in use during an inversion or migration will be
            # overwritten with nonzero values later on.
            for t in d:
                t.data[:] = 0.

            # write traces
            preprocess.writer(d, self.cwd +'/'+ 'traces/adj', filename)


    def check_mesh_properties(self, path=None):
        if not path:
            path = PATH.MODEL_INIT
        if not exists(path):
            raise Exception

        # count slices and grid points
        key = self.parameters[0]
        iproc = 0
        ngll = []
        while True:
            dummy = self.io.read_slice(path, key, iproc)[0]
            ngll += [len(dummy)]
            iproc += 1
            if not exists('%s/proc%06d_%s.bin' % (path, iproc, key)):
                break
        nproc = iproc

        # create coordinate pointers
        coords = Struct()
        for key in ['x', 'y', 'z']:
           coords[key] = partial(self.io.read_slice, self, path, key)

        self._mesh_properties = Struct([
            ['nproc', nproc],
            ['ngll', ngll],
            ['path', path],
            ['coords', coords]])


    def check_source_names(self):
        """ Determines names of sources by applying wildcard rule to user-
            supplied input files
        """
        path = PATH.SPECFEM_DATA
        if not exists(path):
            raise Exception

        # apply wildcard rule
        wildcard = self.source_prefix+'_*'
        globstar = sorted(glob(path +'/'+ wildcard))
        if not globstar:
             print msg.SourceError_SPECFEM % (path, wildcard)
             sys.exit(-1)

        names = []
        for path in globstar:
            names += [basename(path).split('_')[-1]]
        self._source_names = names[:PAR.NTASK]


    def check_solver_parameter_files(self):
        # optional method, can be implemented by subclass
        pass


    ### additional solver attributes

    @property
    def taskid(self):
        # because it is sometimes useful to overload system.taskid
        return system.taskid()

    @property
    def source_name(self):
        # returns name of source currently under consideration
        return self.source_names[self.taskid]

    @property
    def cwd(self):
        # returns working directory currently in use
        return join(PATH.SOLVER, self.source_name)

    @property
    def source_names(self):
       if not hasattr(self, '_source_names'):
           self.check_source_names()
       return self._source_names

    @property
    def mesh_properties(self):
        if not hasattr(self, '_mesh_properties'):
            self.check_mesh_properties()
        return self._mesh_properties

    @property
    def data_filenames(self):
        # required method, must be implemented by subclass
        return NotImplementedError

    @property
    def model_databases(self):
        # required method, must be implemented by subclass
        return NotImplementedError

    @property
    def kernel_databases(self):
        # required method, must be implemented by subclass
        return NotImplementedError

    @property
    def source_prefix(self):
        # required method, must be implemented by subclass
        return NotImplementedError


