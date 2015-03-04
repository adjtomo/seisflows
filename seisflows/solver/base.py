
import subprocess
from glob import glob
from os.path import join

import numpy as np

import seisflows.seistools.specfem3d as solvertools
from seisflows.seistools.shared import getpar, setpar
from seisflows.seistools.io import loadbypar, loadbyproc, savebin, \
    applymap, splitvec, ModelStruct, MinmaxStruct

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import findpath, ParameterObj, ParameterError

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import preprocess


class base(object):
    """ Base class for SPECFEM2D, SPECFEM3D and SPECFEM3D_GLOBE

      eval_func, eval_grad, apply_hess
        These methods deal with evaluation of the misfit function or its
        derivatives.  Together, they provide the primary interface through which
        the solver interacts with other SeisFlows objects.

      forward, adjoint, generate_data, generate_mesh
        These methods allow direct access to low-level SPECFEM components,
        providing another interface through which to interact with the solver.

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

    #  By default, SPECFEM reads velocity models are specified in terms of rho,
    #  vp, and vs. For variable density elastic simulations, include all three
    #  in the 'parameters' list below. For constant density elastic simulations,
    #  remove 'rho' from the list. For variable density acoustic simulations, 
    #  remove 'vs' from the list. For constant desnity acoustic simulations, 
    #  remove 'rho' and 'vs' from the list.

    parameters = []
    parameters += ['vp']
    parameters += ['vs']
    parameters += ['rho']

    # Because density is not well constrained by seismic measurments, it is
    # empirical scaling relations between density and compressional wave 
    # velocity are typically applied.
    density_scaling = None

    #  For seismic inversion, it can be advantageous to express the derivatives
    #  of an objective function in terms of parameters such as bulk c and 
    #  bulk mu. For examples of this approach, see SEISFLOWS-RESEARCH.


    def check(self):
        """ Checks parameters and paths
        """
        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        # check scratch paths
        if 'GLOBAL' not in PATH:
            raise ParameterError(PATH, 'GLOBAL')

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SOLVER' not in PATH:
            if PATH.LOCAL:
                setattr(PATH, 'SOLVER', join(PATH.LOCAL, 'solver'))
            else:
                setattr(PATH, 'SOLVER', join(PATH.GLOBAL, 'solver'))

        # check solver input paths
        if 'SPECFEM_BIN' not in PATH:
            raise ParameterError(PATH, 'SPECFEM_BIN')

        if 'SPECFEM_DATA' not in PATH:
            raise ParameterError(PATH, 'SPECFEM_DATA')


    def setup(self):
        """ Prepares solver for inversion or migration
        """
        unix.rm(self.getpath)

        # As input for an inversion or migration, users can choose between
        # providing data, or providing a target model from which data are
        # generated on the fly. In the former case, a value for PATH.DATA must
        # be supplied; in the latter case, a value for PATH.MODEL_TRUE
        if PATH.DATA:
            self.initialize_solver_directories()
            src = glob(PATH.DATA +'/'+ self.getname +'/'+ '*')
            dst = 'traces/obs/'
            unix.cp(src, dst)

        else:
            self.generate_data(
                model_path=PATH.MODEL_TRUE,
                model_name='model_true',
                model_type='gll')

        # prepare model
        self.generate_mesh(
            model_path=PATH.MODEL_INIT,
            model_name='model_init',
            model_type='gll')

        self.initialize_adjoint_traces()
        self.initialize_io_machinery()


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
        unix.mv(self.data_wildcard, 'traces/syn')
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
        """ Computes action of Hessian on a given model vector.
        """
        unix.cd(self.getpath)
        unix.mkdir('traces/lcg')

        self.import_model(path)
        self.forward()
        unix.mv(self.data_wildcard, 'traces/lcg')
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

    def load(self, path, parameters=None, mapping=None, prefix='', suffix='', verbose=False):
        """ reads SPECFEM model

          Models are stored in Fortran binary format and separated into multiple
          files according to material parameter and processor rank.

          Optionally, 'mapping' can be used to convert on the fly from one set
          of material parameters to another.
        """
        if not parameters:
            parameters = self.parameters

        model = ModelStruct(parameters, mapping)
        minmax = MinmaxStruct(parameters, mapping)

        for iproc in range(PAR.NPROC):
            keys, vals = loadbypar(path, parameters, iproc, prefix, suffix)

            # keep track of min, max
            minmax.update(keys, vals)

            # apply optional mapping
            if mapping:
                keys, vals = applymap(vals, mapping)

            # append latest values
            for key, val in zip(keys, vals):
                model[key] += [val]

        if verbose:
            minmax.write(path, logpath=PATH.SUBMIT)

        return model


    def save(self, path, model, prefix='', suffix=''):
        """ writes SPECFEM model

            The following code writes SPECFEM acoustic and elastic models.
            For code that writes transversely isotropic models, see 
            solver.specfem3d_globe.

            There is a large tradeoff here between being simple and being 
            flexible.  In this case we opt for a simple hardwired approach. For
            a much more flexible approach, see seisflows-research.
        """
        unix.mkdir(path)

        for iproc in range(PAR.NPROC):
            if 'vp' in self.parameters:
                savebin(model['vp'][iproc], path, iproc, 'vp')
            else:
                self.copybin(path, iproc, 'vp')

            if 'vs' in self.parameters:
                savebin(model['vs'][iproc], path, iproc, 'vs')
            else:
                self.copybin(path, iproc, 'vs')

            if 'rho' in self.parameters:
                savebin(model['rho'][iproc], path, iproc, 'rho')
            elif self.density_scaling:
                raise NotImplementedError
            else:
                self.copybin(path, iproc, 'rho')


    def merge(self, model):
        """ Converts model from dictionary to vector representation
        """
        v = np.array([])
        for key in self.parameters:
            for iproc in range(PAR.NPROC):
                v = np.append(v, model[key][iproc])
        return v


    def split(self, v):
        """ Converts model from vector to dictionary representation
        """
        nproc = PAR.NPROC
        ndim = len(self.parameters)
        npts = len(v)/(nproc*ndim)

        model = {}
        for idim, key in enumerate(self.parameters):
            model[key] = splitvec(v, nproc, npts, idim)
        return model



    ### postprocessing utilities

    def combine(self, path=''):
        """ Sums individual source contributions. Wrapper over xcombine_sem
            utility.
        """
        unix.cd(self.getpath)

        with open('kernel_paths', 'w') as f:
            f.writelines([join(path, dir)+'\n' for dir in unix.ls(path)])

        unix.mkdir(path +'/'+ 'sum')
        for name in self.parameters:
            self.mpirun(
                PATH.SPECFEM_BIN +'/'+ 'xcombine_sem '
                + name + '_kernel' + ' '
                + 'kernel_paths' + ' '
                + path +'/'+ 'sum')


    def smooth(self, path='', span=0.):
        """ Smooths kernels by convolving them with a Gaussian.  Wrapper over 
            xsmooth_sem utility.
        """
        assert (exists(path))

        # apply smoothing operator
        unix.cd(self.getpath)
        for name in self.parameters:
            print ' smoothing', name
            self.mpirun(
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


    def clip(self, path='', minval=-999999., maxval=999999., thresh='dummy'):
        """ Clips kernels by convolving them with a Gaussian.  Wrapper over 
            xclip_sem utility.
        """
        assert (exists(path))

        # apply smoothing operator
        unix.cd(self.getpath)
        for name in self.parameters:
            self.mpirun(
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

        if system.getnode()==0:
            self.save(dst, self.load(src, verbose=True))
        else:
            self.save(dst, self.load(src))

    def import_traces(self, path):
        src = glob(join(path, 'traces', self.getname, '*'))
        dst = join(self.getpath, 'traces/obs')
        unix.cp(src, dst)

    def export_model(self, path):
        if system.getnode() == 0:
            src = glob(join(self.model_databases, '*.bin'))
            dst = path
            unix.mkdir(dst)
            unix.cp(src, dst)

    def export_kernels(self, path):
        # work around unfortunate SPECFEM conventions
        try:
            files = glob(self.model_databases +'/'+ '*alpha*_kernel.bin')
            unix.rename('alpha', 'vp', files)

            files = glob(self.model_databases +'/'+ '*beta*_kernel.bin')
            unix.rename('beta', 'vs', files)
        except:
            pass

        # export kernels
        unix.mkdir_gpfs(join(path, 'kernels'))
        unix.mkdir(join(path, 'kernels', self.getname))
        src = join(glob(self.model_databases +'/'+ '*kernel.bin'))
        dst = join(path, 'kernels', self.getname)
        unix.mv(src, dst)

    def export_residuals(self, path):
        unix.mkdir_gpfs(join(path, 'residuals'))
        src = join(self.getpath, 'residuals')
        dst = join(path, 'residuals', self.getname)
        unix.mv(src, dst)

    def export_traces(self, path, prefix='traces/obs'):
        unix.mkdir_gpfs(join(path, 'traces'))
        src = join(self.getpath, prefix)
        dst = join(path, 'traces', self.getname)
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

        unix.mkdir('traces/obs')
        unix.mkdir('traces/syn')
        unix.mkdir('traces/adj')
        unix.mkdir(self.model_databases)

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


    def initialize_adjoint_traces(self):
        """ Adjoint traces are initialized by writing zeros for all components.
            Components actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        _, h = preprocess.load('traces/obs')
        zeros = np.zeros((h.nt, h.nr))
        for channel in ['x', 'y', 'z']:
            preprocess.writer(zeros, h, channel=channel, prefix='traces/adj/')


    def initialize_io_machinery(self):
        """ Writes mesh files expected by input/output methods
        """
        if system.getnode() == 0:
            if 'OPTIMIZE' in PATH:
                model = self.load(PATH.MODEL_INIT)
                if not exists(PATH.OPTIMIZE +'/'+ 'm_new'):
                    savenpy(PATH.OPTIMIZE +'/'+ 'm_new', self.merge(model))


    ### miscellaneous

    def mpirun(self, script, output='/dev/null'):
        """ Wrapper for mpirun
        """
        with open(output,'w') as f:
            subprocess.call(
                system.mpiargs() + script,
                shell=True,
                stdout=f)


    def copybin(self, path, iproc, par):
        """ Copies SPECFEM database file
        """
        src = PATH.OUTPUT +'/'+ 'model_init'
        dst = path

        filename = 'proc%06d_%s.bin' % (iproc, par)
        unix.cp(join(src, filename), join(dst, filename))


    @property
    def getname(self):
        """ Returns name of source currently under consideration
        """
        isrc = system.getnode()
        if not hasattr(self, 'sources'):
            paths = glob(PATH.SPECFEM_DATA +'/'+ self.source_prefix+'_*')
            self.sources = []
            for path in paths:
                self.sources += [unix.basename(path).split('_')[-1]]
            self.sources.sort()
        return self.sources[isrc]

    @property
    def getpath(self):
        """ Returns working directory corresponding to current source
        """
        return join(PATH.SOLVER, self.getname)


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

