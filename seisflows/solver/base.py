
import subprocess
import sys
import numpy as np

from functools import partial
from glob import glob
from os.path import basename, join
from seisflows.config import ParameterError, custom_import
from seisflows.plugins.io import sem
from seisflows.tools import msg, unix
from seisflows.tools.shared import ModelDict
from seisflows.tools.tools import Struct, diff, exists, call_solver



PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']



class base(object):
    """ Base class for SPECFEM2D, SPECFEM3D and SPECFEM3D_GLOBE

      eval_func, eval_grad, apply_hess
        These methods deal with evaluation of the misfit function or its
        derivatives.  Together, they provide the primary interface through which
        SeisFlows interacts with SPECFEM.

      forward, adjoint
        These methods allow direct access to low-level SPECFEM components,
        providing another interface through which to interact with the solver.

     generate_data, generate_model
        One time operations performed at beginning of an inversion or 
        migration.

     initialize_solver_directories, initialize_adjoint_traces
        SPECFEM requires a particular directory structure in which to run and
        particular file formats for models, data, and parameter files. These
        methods help put in place all these prerequisites.

      load, save
        For reading and writing SPECFEM models and kernels. On the disk,
        models and kernels are stored as binary files, and in memory, as
        dictionaries with different keys corresponding to different material
        parameters.

      split, merge
        Within the solver routines, it is natural to store models as 
        dictionaries. Within the optimization routines, it is natural to store
        models as vectors. Two methods, 'split' and 'merge', are used to convert 
        back and forth between these two representations.

      combine, smooth
        Utilities for combining and smoothing kernels.
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
        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        # check scratch paths
        if 'SCRATCH' not in PATH:
            raise ParameterError(PATH, 'SCRATCH')

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SOLVER' not in PATH:
            if PATH.LOCAL:
                setattr(PATH, 'SOLVER', join(PATH.LOCAL, 'solver'))
            else:
                setattr(PATH, 'SOLVER', join(PATH.SCRATCH, 'solver'))

        # check solver input paths
        if 'SPECFEM_BIN' not in PATH:
            raise ParameterError(PATH, 'SPECFEM_BIN')

        if 'SPECFEM_DATA' not in PATH:
            raise ParameterError(PATH, 'SPECFEM_DATA')

        # assertions
        assert self.parameters != []


    def setup(self):
        """ Prepares solver for inversion or migration
        """
        # clean up for new inversion
        unix.rm(self.getpath)

        # As input for an inversion or migration, users can choose between
        # providing data, or providing a target model from which data are
        # generated on the fly. In the former case, a value for PATH.DATA must
        # be supplied; in the latter case, a value for PATH.MODEL_TRUE must be
        # provided

        if PATH.DATA:
            # copy user supplied data
            self.initialize_solver_directories()

            src = glob(PATH.DATA +'/'+ self.getname +'/'+ '*')
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
        unix.cd(self.getpath)
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

    def eval_func(self, path='', export_traces=False):
        """ Evaluates misfit function by carrying out forward simulation and
            comparing observations and synthetics.
        """
        unix.cd(self.getpath)
        self.import_model(path)

        self.forward()
        preprocess.prepare_eval_grad(self.getpath)
        self.export_residuals(path)


    def eval_grad(self, path='', export_traces=False):
        """ Evaluates gradient by carrying out adjoint simulation. Adjoint traces
            must be in place beforehand.
        """
        unix.cd(self.getpath)

        self.adjoint()

        self.export_kernels(path)

        if export_traces:
            self.export_traces(path+'/'+'traces/syn', prefix='traces/syn')
            self.export_traces(path+'/'+'traces/adj', prefix='traces/adj')


    def apply_hess(self, path=''):
        """ Computes action of Hessian on a given model vector. A gradient 
          evaluation must have already been carried out.
        """
        unix.cd(self.getpath)
        unix.mkdir('traces/lcg')

        self.import_model(path)
        self.forward('traces/lcg')
        preprocess.prepare_apply_hess(self.getpath)

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

    def load(self, path, parameters=[], prefix='', suffix=''):
        """ reads SPECFEM model or kernels
        """
        dict = ModelDict()
        for iproc in range(self.mesh_properties.nproc):
            for key in parameters or self.parameters:
                dict[key] += sem.read(path, prefix+key+suffix, iproc)
        return dict


    def save(self, dict, path, parameters=['vp','vs','rho'], prefix='', suffix=''):
        """ writes SPECFEM model or kernels
        """
        unix.mkdir(path)

        # fill in missing keys
        missing_keys = diff(parameters, dict.keys())
        for iproc in range(self.mesh_properties.nproc):
            for key in missing_keys:
                dict[key] += sem.read(PATH.MODEL_INIT, prefix+key+suffix, iproc)

        # save values
        for iproc in range(self.mesh_properties.nproc):
            for key in parameters:
                sem.write(dict[key][iproc], path, prefix+key+suffix, iproc)


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
        model = ModelDict()
        for idim, key in enumerate(parameters or self.parameters):
            model[key] = []
            for iproc in range(nproc):
                imin = sum(ngll)*idim + sum(ngll[:iproc])
                imax = sum(ngll)*idim + sum(ngll[:iproc+1])
                model[key] += [m[imin:imax]]
        return model



    ### postprocessing wrappers

    def combine(self, path='', parameters=[]):
        """ Sums individual source contributions. Wrapper over xcombine_sem
            utility.
        """
        unix.cd(self.getpath)

        with open('kernel_paths', 'w') as file:
            file.writelines([join(path, dir)+'\n' for dir in self.source_names])

        unix.mkdir(path +'/'+ 'sum')
        for name in parameters or self.parameters:
            call_solver(
                system.mpiexec(),
                PATH.SPECFEM_BIN +'/'+ 'xcombine_sem '
                + name + '_kernel' + ' '
                + 'kernel_paths' + ' '
                + path +'/'+ 'sum')


    def smooth(self, path='', parameters=[], span=0.):
        """ Smooths kernels by convolving them with a Gaussian.  Wrapper over 
            xsmooth_sem utility.
        """
        assert exists(path)
        assert len(parameters) > 0

        # apply smoothing operator
        unix.cd(self.getpath)
        for name in parameters or self.parameters:
            print ' smoothing', name
            call_solver(
                system.mpiexec(),
                PATH.SPECFEM_BIN +'/'+ 'xsmooth_sem '
                + str(span) + ' '
                + str(span) + ' '
                + name + '_kernel' + ' '
                + path + '/ '
                + path + '/ ',
                output='/dev/null')

        print ''

        # move input files
        src = path
        dst = path + '_nosmooth'
        unix.mkdir(dst)
        for name in parameters or self.parameters:
            unix.mv(glob(src+'/*'+name+'_kernel.bin'), dst)

        # rename output files
        unix.rename('_smooth', '', glob(src+'/*'))


    def clip(self, path='', parameters=[], minval=-np.inf, maxval=np.inf):
        """ Clips kernels by convolving them with a Gaussian.  Wrapper over 
            xclip_sem utility.
        """
        assert exists(path)
        assert len(parameters) > 0

        unix.cd(self.getpath)
        for name in parameters or self.parameters:
            call_solver(
                system.mpiexec,
                PATH.SPECFEM_BIN +'/'+ 'xclip_sem '
                + str(minval) + ' '
                + str(maxval) + ' '
                + name + '_kernel' + ' '
                + path + '/ '
                + path + '/ ')

        # move input files
        src = path
        dst = path + '_noclip'
        unix.mkdir(dst)
        for name in parameters or self.parameters:
            unix.mv(glob(src+'/*'+name+'.bin'), dst)

        # rename output files
        unix.rename('_clip', '', glob(src+'/*'))


    ### file transfer utilities

    def import_model(self, path):
        model = self.load(path+'/'+'model')
        self.save(model, self.model_databases)

    def import_traces(self, path):
        src = glob(join(path, 'traces', self.getname, '*'))
        dst = join(self.getpath, 'traces/obs')
        unix.cp(src, dst)

    def export_model(self, path, solver_parameters=['rho', 'vp', 'vs']):
        if self.getnode == 0:
            unix.mkdir(path)
            for key in solver_parameters:
                files = glob(join(self.model_databases, '*'+key+'.bin'))
                unix.cp(files, path)

    def export_kernels(self, path):
        unix.cd(self.kernel_databases)

        # work around conflicting name conventions
        self.rename_kernels()

        src = glob('*_kernel.bin')
        dst = join(path, 'kernels', self.getname)
        unix.mkdir(dst)
        unix.mv(src, dst)

    def export_residuals(self, path):
        unix.mkdir(join(path, 'residuals'))

        src = join(self.getpath, 'residuals')
        dst = join(path, 'residuals', self.getname)
        unix.mv(src, dst)

    def export_traces(self, path, prefix='traces/obs'):
        unix.mkdir(join(path))

        src = join(self.getpath, prefix)
        dst = join(path, self.getname)
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
        unix.mkdir(self.getpath)
        unix.cd(self.getpath)

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

        src = 'DATA/' + self.source_prefix +'_'+ self.getname
        dst = 'DATA/' + self.source_prefix
        unix.cp(src, dst)

        self.check_solver_parameter_files()


    def initialize_adjoint_traces(self):
        """ Adjoint traces are initialized by writing zeros for all components.
            Components actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        for filename in self.data_filenames:
            # read traces
            d = preprocess.reader(self.getpath +'/'+ 'traces/obs', filename)

            # replace data with zeros
            for t in d:
                t.data[:] = 0.

            # write traces
            preprocess.writer(d, self.getpath +'/'+ 'traces/adj', filename)


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
            dummy = sem.read(path, key, iproc)[0]
            ngll += [len(dummy)]
            iproc += 1
            if not exists('%s/proc%06d_%s.bin' % (path, iproc, key)):
                break
        nproc = iproc

        # create coordinate pointers
        coords = Struct()
        for key in ['x', 'y', 'z']:
           coords[key] = partial(sem.read, path, key)

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

        # apply wildcard
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
    def getnode(self):
        # because it is sometimes useful to overload system.getnode
        return system.getnode()

    @property
    def getname(self):
        # returns name of source currently under consideration
        return self.source_names[self.getnode]

    @property
    def getpath(self):
        # returns working directory currently in use
        return join(PATH.SOLVER, self.getname)

    @property
    def mesh_properties(self):
        if not hasattr(self, '_mesh_properties'):
            self.check_mesh_properties()
        return self._mesh_properties

    @property
    def source_names(self):
       if not hasattr(self, '_source_names'):
           self.check_source_names()
       return self._source_names

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

