
import subprocess
import sys
from glob import glob
from os.path import basename, join

import numpy as np

import seisflows.seistools.specfem3d as solvertools
from seisflows.seistools.shared import getpar, setpar, Model, Minmax
from seisflows.seistools.io import loadbypar, copybin, loadbin, savebin, splitvec

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.code import Struct, exists, mpicall
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import preprocess


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
        density_scaling = None
        parameters += ['rho']
    elif PAR.DENSITY == 'Constant':
        density_scaling = None


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

            src = glob(PATH.DATA +'/'+ basename(self.getpath) +'/'+ '*')
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
        pass
        #unix.rm('OUTPUT_FILES')
        #unix.mkdir('OUTPUT_FILES')


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
        unix.mv(self.getdata, 'traces/syn')
        preprocess.prepare_eval_grad(self.getpath)

        self.export_residuals(path)
        if export_traces:
            self.export_traces(path, prefix='traces/syn')


    def eval_grad(self, path='', export_traces=False):
        """ Evaluates gradient by carrying out adjoint simulation. Adjoint traces
            must be in place beforehand.
        """
        unix.cd(self.getpath)

        self.adjoint()

        self.export_kernels(path)
        if export_traces:
            self.export_traces(path, prefix='traces/adj')


    def apply_hess(self, path=''):
        """ Computes action of Hessian on a given model vector. A gradient 
          evaluation must have already been carried out beforehand.
        """
        unix.cd(self.getpath)
        unix.mkdir('traces/lcg')

        self.import_model(path)
        self.forward()
        unix.mv(self.channels, 'traces/lcg')
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

    def load(self, path, prefix='', suffix='', verbose=False):
        """ reads SPECFEM model or kernel

          Models are stored in Fortran binary format and separated into multiple
          files according to material parameter and processor rank.
        """
        minmax = Minmax(self.parameters)
        model = Model(self.parameters)

        for iproc in range(self.mesh.nproc):
            # read database files
            keys, vals = loadbypar(path, self.parameters, iproc, prefix, suffix)
            for key, val in zip(keys, vals):
                model[key] += [val]

            minmax.update(keys, vals)

        if verbose:
            minmax.write(path, logpath=PATH.SUBMIT)

        return model


    def save(self, path, model, prefix='', suffix='', solver_parameters=['vp', 'vs']):
        """ writes SPECFEM model

            The following method writes acoustic and elastic models. For an
            example of how to write transversely isotropic models, see 
            solver.specfem3d_globe.

            There is a tradeoff here between being simple and being 
            flexible.  In this case we opt for a simple hardwired approach. For
            a more flexible, more complex approach see SEISFLOWS-RESEARCH.
        """
        unix.mkdir(path)

        if 'kernel' in suffix:
            for iproc in range(self.mesh.nproc):
                for key in solver_parameters:
                    if key in self.parameters:
                        savebin(model[key][iproc], path, iproc, prefix+key+suffix)
                if 'rho' in self.parameters:
                    savebin(model['rho'][iproc], path, iproc, prefix+'rho'+suffix)

        else:
            for iproc in range(self.mesh.nproc):
                for key in solver_parameters:
                    if key in self.parameters:
                        savebin(model[key][iproc], path, iproc, prefix+key+suffix)
                    else:
                        src = PATH.OUTPUT +'/'+ 'model_init'
                        dst = path
                        copybin(src, dst, iproc, prefix+key+suffix)

                if 'rho' in self.parameters:
                    savebin(model['rho'][iproc], path, iproc, prefix+'rho'+suffix)
                else:
                    src = PATH.OUTPUT +'/'+ 'model_init'
                    dst = path
                    copybin(src, dst, iproc, 'rho')



    def merge(self, model):
        """ Converts model from dictionary to vector representation
        """
        v = np.array([])
        for key in self.parameters:
            for iproc in range(self.mesh.nproc):
                v = np.append(v, model[key][iproc])
        return v


    def split(self, v):
        """ Converts model from vector to dictionary representation
        """
        model = {}
        for idim, key in enumerate(self.parameters):
            model[key] = splitvec(v, self.mesh.nproc, self.mesh.ngll, idim)
        return model



    ### postprocessing utilities

    def combine(self, path='', parameters=[]):
        """ Sums individual source contributions. Wrapper over xcombine_sem
            utility.
        """
        unix.cd(self.getpath)

        names = self.check_source_names()
        with open('kernel_paths', 'w') as f:
            f.writelines([join(path, dir)+'\n' for dir in names])

        unix.mkdir(path +'/'+ 'sum')
        for name in parameters:
            mpicall(
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
        for name in parameters:
            print ' smoothing', name
            mpicall(
                system.mpiexec(),
                PATH.SPECFEM_BIN +'/'+ 'xsmooth_sem '
                + str(span) + ' '
                + str(span) + ' '
                + name + '_kernel' + ' '
                + path + '/ '
                + path + '/ ',
                output=self.getpath+'/'+'OUTPUT_FILES/output_smooth_sem.txt')

        print ''

        # move input files
        src = path
        dst = path + '_nosmooth'
        unix.mkdir(dst)
        for name in self.parameters:
            unix.mv(glob(src+'/*'+name+'.bin'), dst)

        # rename output files
        unix.rename('_smooth', '', glob(src+'/*'))


    def clip(self, path='', parameters=[], minval=-np.inf, maxval=np.inf):
        """ Clips kernels by convolving them with a Gaussian.  Wrapper over 
            xclip_sem utility.
        """
        assert exists(path)
        assert len(parameters) > 0

        unix.cd(self.getpath)
        for name in self.parameters:
            mpicall(
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
        for name in self.parameters:
            unix.mv(glob(src+'/*'+name+'.bin'), dst)

        # rename output files
        unix.rename('_clip', '', glob(src+'/*'))


    ### file transfer utilities

    def import_model(self, path):
        src = join(path, 'model')
        dst = self.model_databases

        if self.getnode==0:
            self.save(dst, self.load(src, verbose=True))
        else:
            self.save(dst, self.load(src))

    def import_traces(self, path):
        src = glob(join(path, 'traces', basename(self.getpath), '*'))
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

        # hack to deal with problems on parallel filesystem
        unix.mkdir(join(path, 'kernels'), noexit=True)

        unix.mkdir(join(path, 'kernels', basename(self.getpath)))
        src = join(glob('*_kernel.bin'))
        dst = join(path, 'kernels', basename(self.getpath))
        unix.mv(src, dst)

    def export_residuals(self, path):
        # hack deals with problems on parallel filesystem
        unix.mkdir(join(path, 'residuals'), noexit=True)

        src = join(self.getpath, 'residuals')
        dst = join(path, 'residuals', basename(self.getpath))
        unix.mv(src, dst)

    def export_traces(self, path, prefix='traces/obs'):
        # hack deals with problems on parallel filesystem
        unix.mkdir(join(path, 'traces'), noexit=True)

        src = join(self.getpath, prefix)
        dst = join(path, 'traces', basename(self.getpath))
        unix.cp(src, dst)


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

        src = 'DATA/' + self.source_prefix +'_'+ basename(self.getpath)
        dst = 'DATA/' + self.source_prefix
        unix.cp(src, dst)

        self.check_solver_parameter_files()


    def initialize_adjoint_traces(self):
        """ Adjoint traces are initialized by writing zeros for all components.
            Components actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        path = self.getpath+'/'+'traces/adj/'
        for channel in ['x', 'y', 'z']:
            preprocess.write_zero_traces(path, channel)


    def check_mesh_properties(self, path=None, parameters=None):
        if not hasattr(self, '_mesh_properties'):
            if not path:
                path = PATH.MODEL_INIT

            if not parameters:
                parameters = self.parameters

            nproc = 0
            ngll = []
            while True:
                dummy = loadbin(path, nproc, parameters[0])
                ngll += [len(dummy)]
                nproc += 1
                if not exists('%s/proc%06d_%s.bin' % (path, nproc, parameters[0])):
                    break

            self._mesh_properties = Struct([
                ['nproc', nproc],
                ['ngll', ngll]])

        return self._mesh_properties


    def check_source_names(self):
        """ Checks names of sources
        """
        if not hasattr(self, '_source_names'):
            path = PATH.SPECFEM_DATA
            wildcard = self.source_prefix+'_*'
            globstar = sorted(glob(path +'/'+ wildcard))
            if not globstar:
                 print msg.SourceError_SPECFEM % (path, wildcard)
                 sys.exit(-1)
            names = []
            for path in globstar:
                names += [basename(path).split('_')[-1]]
            self._source_names = names[:PAR.NTASK]

        return self._source_names


    def check_solver_parameter_files(self):
        # must be implemented by subclass
        pass


    ### solver utility functions

    @property
    def getnode(self):
        # it is sometimes useful to overload system.getnode via this method
        return system.getnode()

    @property
    def getname(self):
        # returns name of source currently under consideration
        return self.check_source_names()[self.getnode]

    @property
    def getpath(self):
        # returns solver working directory currently in use
        return join(PATH.SOLVER, self.getname)

    @property
    def getdata(self):
        # returns data filenames for current source
        wildcard = join(self.getpath, 'OUTPUT_FILES', self.data_wildcard)
        return glob(wildcard)


    ### solver information

    @property
    def mesh(self):
        return self.check_mesh_properties()

    @property
    def channels(self):
        wildcard = join(self.getpath, 'traces/obs', self.data_wildcard)
        return [basename(filename) for filename in glob(wildcard)]


    # solver paths and filenames

    @property
    def data_wildcard(self):
        # must be implemented by subclass
        return NotImplementedError

    @property
    def model_databases(self):
        # must be implemented by subclass
        return NotImplementedError

    @property
    def source_prefix(self):
        # must be implemented by subclass
        return NotImplementedError

