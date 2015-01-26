
import subprocess
from os.path import join
from glob import glob

import numpy as np

import seisflows.seistools.specfem2d as solvertools

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists, setdiff
from seisflows.tools.config import findpath, ParameterObj
from seisflows.tools.io import loadbin, savebin

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import preprocess


class specfem2d(object):
    """ Python interface for SPECFEM2D

      eval_func, eval_grad, apply_hess
        These methods deal with evaluation of the misfit function or its
        derivatives and provide the primary interface between the solver and
        other workflow components.

      forward, adjoint
        These methods allow direct access to individual SPECFEM2D components.
        Together, they provide a secondary interface users can employ for
        specialized tasks not covered by high level methods.

      prepare_solver, prepare_data, prepare_model
        SPECFEM2D requires a particular directory structure in which to run and
        particular file formats for models, data, and parameter files. These
        methods help put in place all these prerequisites.

      load, save
        For reading and writing SPECFEM2D models and kernels. On the disk,
        models and kernels are stored as binary files, and in memory, as
        dictionaries with different keys corresponding to different material
        parameters.

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

    kernel_map = {
        'rho': 'rho_kernel',
        'vp': 'alpha_kernel',
        'vs': 'beta_kernel'}

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # check time stepping parameters
        if 'NT' not in PAR:
            raise Exception

        if 'DT' not in PAR:
            raise Exception

        if 'F0' not in PAR:
            raise Exception

        # check paths
        if 'GLOBAL' not in PATH:
            raise Exception

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SOLVER' not in PATH:
            if PATH.LOCAL:
                setattr(PATH, 'SOLVER', join(PATH.LOCAL, 'solver'))
            else:
                setattr(PATH, 'SOLVER', join(PATH.GLOBAL, 'solver'))


    def setup(self):
        """ Prepares solver for inversion or migration

          As input for an inversion or migration, users can choose between
          supplying data or providing a target model from which data are
          generated on the fly. In both cases, all necessary SPECFEM2D input
          files must be provided.
        """
        unix.rm(self.getpath)

        # prepare data
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

    def generate_data(self, **model_kwargs):
        """ Generates data
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.getpath)
        solvertools.setpar('SIMULATION_TYPE', '1')
        solvertools.setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D')

        unix.mv(self.wildcard, 'traces/obs')
        self.export_traces(PATH.OUTPUT, 'traces/obs')


    def generate_mesh(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(model_type)
        assert (exists(model_path))

        self.initialize_solver_directories()
        unix.cp(model_path, 'DATA/model_velocity.dat_input')
        self.export_model(PATH.OUTPUT +'/'+ model_name)


    ### high-level solver interface

    def eval_func(self, path='', export_traces=False):
        """ Evaluates misfit function by carrying out forward simulation and
            making measurements on observations and synthetics.
        """
        unix.cd(self.getpath)
        self.import_model(path)

        self.forward()
        unix.mv(self.wildcard, 'traces/syn')
        preprocess.prepare_eval_grad(self.getpath)

        self.export_residuals(path)
        if export_traces:
            self.export_traces(path, prefix='traces/syn')


    def eval_grad(self, path='', export_traces=False):
        """ Evaluates gradient by carrying out adjoint simulation. Adjoint traces
            must be in place prior to calling this method.
        """
        unix.cd(self.getpath)

        self.adjoint()

        self.export_kernels(path)
        if export_traces:
            self.export_traces(path, prefix='traces/syn')


    def apply_hess(self, path='', hessian='Newton'):
        """ Evaluates action of Hessian on a given model vector.
        """
        unix.cd(self.getpath)
        self.imprt(path, 'model')

        self.forward()
        unix.mv(self.wildcard, 'traces/lcg')
        preprocess.prepare_apply_hess(self.getpath)
        self.adjoint()

        self.export_kernels(path)



    ### low-level solver interface

    def forward(self):
        """ Calls SPECFEM2D forward solver
        """
        solvertools.setpar('SIMULATION_TYPE', '1')
        solvertools.setpar('SAVE_FORWARD', '.true.')

        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D')
        unix.mv(self.wildcard, 'traces/syn')


    def adjoint(self):
        """ Calls SPECFEM2D adjoint solver
        """
        solvertools.setpar('SIMULATION_TYPE', '3')
        solvertools.setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')

        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D')



    ### model input/output

    def load(self, filename, type='', verbose=False):
        """Reads SPECFEM2D kernel or model

           Models and kernels are read from 5 or 6 column text files whose
           format is described in the SPECFEM2D user manual. Once read, a model
           or kernel is stored in a dictionary containing mesh coordinates and
           corresponding material parameter values.
        """
        # read text file
        M = np.loadtxt(filename)
        nrow = M.shape[0]
        ncol = M.shape[1]

        if ncol == 5:
            ioff = 0
        elif ncol == 6:
            ioff = 1
        else:
            raise ValueError("Wrong number of columns.")

        # fill in dictionary
        parts = {}
        for key in ['x', 'z', 'rho', 'vp', 'vs']:
            parts[key] = [M[:,ioff]]
            ioff += 1
        return parts

    def save(self, filename, parts, type='model'):
        """writes SPECFEM2D kernel or model"""
        # allocate array
        if type == 'model':
            nrow = len(parts[parts.keys().pop()][0])
            ncol = 6
            ioff = 1
            M = np.zeros((nrow, ncol))
        elif type == 'kernel':
            nrow = len(parts[parts.keys().pop()][0])
            ncol = 5
            ioff = 0
            M = np.zeros((nrow, ncol))
        else:
            raise ValueError

        # fill in array
        for icol, key in enumerate(['x', 'z', 'rho', 'vp', 'vs']):
            M[:,icol+ioff] = parts[key][0]

        # write array
        np.savetxt(filename, M, '%10.4e')



    ### vector/dictionary conversion

    def merge(self, parts):
        """ merges dictionary into vector
        """
        v = np.array([])
        for key in self.inversion_parameters:
            for iproc in range(PAR.NPROC):
                v = np.append(v, parts[key][iproc])
        return v


    def split(self, v):
        """ splits vector into dictionary
        """
        parts = {}
        nrow = len(v)/(PAR.NPROC*len(self.inversion_parameters))
        j = 0
        for key in ['x', 'z', 'rho', 'vp', 'vs']:
            parts[key] = []
            if key in self.inversion_parameters:
                for i in range(PAR.NPROC):
                    imin = nrow*PAR.NPROC*j + nrow*i
                    imax = nrow*PAR.NPROC*j + nrow*(i + 1)
                    i += 1
                    parts[key].append(v[imin:imax])
                j += 1
            else:
                for i in range(PAR.NPROC):
                    proc = '%06d' % i
                    part = np.load(PATH.GLOBAL +'/'+ 'mesh' +'/'+ key +'/'+ proc)
                    parts[key].append(part)
        return parts



    ### postprocessing utilities

    def combine(self, path=''):
        """combines SPECFEM2D kernels"""
        subprocess.call(
            [self.getpath +'/'+ 'bin/xsmooth_sem'] +
            [str(len(unix.ls(path)))] +
            [path])

    def smooth(self, path='', tag='gradient', span=0.):
        """smooths SPECFEM2D kernels by convolving them with a Gaussian"""
        from seisflows.tools.array import meshsmooth

        parts = self.load(path +'/'+ tag)
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
            parts[key] = [meshsmooth(x, z, parts[key][0], span, nx, nz)]
        unix.mv(path +'/'+ tag, path +'/'+ '_nosmooth')
        self.save(path +'/'+ tag, parts)


    def clip(self, path='', tag='gradient', thresh=1.):
        """clips SPECFEM2D kernels"""
        parts = self.load(path +'/'+ tag)
        if thresh >= 1.:
            return parts

        for key in self.inversion_parameters:
            # scale to [-1,1]
            minval = parts[key][0].min()
            maxval = parts[key][0].max()
            np.clip(parts[key][0], thresh*minval, thresh*maxval, out=parts[key][0])
        unix.mv(path +'/'+ tag, path +'/'+ '_noclip')
        self.save(path +'/'+ tag, parts)


    ### file transfer utilities

    def import_model(self, path):
        src = join(path +'/'+ 'model')
        dst = join(self.getpath, 'DATA/model_velocity.dat_input')
        unix.cp(src, dst)

    def import_traces(self, path):
        src = glob(join(path, 'traces', self.getname, '*'))
        dst = join(self.getpath, 'traces/obs')
        unix.cp(src, dst)

    def export_model(self, path):
        if system.getnode() == 0:
            src = join(self.getpath, 'DATA/model_velocity.dat_input')
            dst = path
            unix.cp(src, dst)

    def export_kernels(self, path):
        unix.mkdir_gpfs(join(path, 'kernels'))
        src = join(self.getpath, 'OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat')
        dst = join(path, 'kernels', '%06d' % system.getnode())
        unix.cp(src, dst)

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
        """ Creates directory structure expected by SPECFEM2D, copies 
          executables, and prepares input files. Executables must be supplied 
          by user as there is currently no mechanism to automatically compile 
          from source.
        """
        unix.mkdir(self.getpath)
        unix.cd(self.getpath)

        # create directory structure
        unix.mkdir('bin')
        unix.mkdir('DATA')

        unix.mkdir('traces/obs')
        unix.mkdir('traces/syn')
        unix.mkdir('traces/adj')

        unix.mkdir(self.databases)

        # copy exectuables
        src = glob(PATH.SOLVER_BINARIES +'/'+ '*')
        dst = 'bin/'
        unix.cp(src, dst)

        # copy input files
        src = glob(PATH.SOLVER_FILES +'/'+ '*')
        dst = 'DATA/'
        unix.cp(src, dst)

        src = 'DATA/SOURCE_' + self.getname
        dst = 'DATA/SOURCE'
        unix.cp(src, dst)

        solvertools.setpar('f0', PAR.F0, 'DATA/SOURCE')


    def initialize_adjoint_traces(self):
        """ Adjoint traces must be initialized by writing zeros for all 
          components. This is because when reading traces at the start of an
          adjoint simulation, SPECFEM2D expects that all components exist.
          Components actually in use during an inversion or migration will
          be overwritten with nonzero values later on.
        """
        _, h = preprocess.load('traces/obs')
        zeros = np.zeros((h.nt, h.nr))
        for channel in ['x', 'y', 'z']:
            preprocess.writer(zeros, h, channel=channel, prefix='traces/adj')


    def initialize_io_machinery(self):
        """ Writes mesh files expected by input/output methods
        """
        if system.getnode() == 0:

            model_set = set(self.model_parameters)
            inversion_set = set(self.inversion_parameters)

            parts = self.load(PATH.MODEL_INIT)
            try:
                path = PATH.GLOBAL +'/'+ 'mesh'
            except:
                raise Exception
            if not exists(path):
                for key in list(setdiff(model_set, inversion_set)) + ['x', 'z']:
                    unix.mkdir(path +'/'+ key)
                    for proc in range(PAR.NPROC):
                        with open(path +'/'+ key +'/'+ '%06d' % proc, 'w') as file:
                            np.save(file, parts[key][proc])

            try:
                path = PATH.OPTIMIZE +'/'+ 'm_new'
            except:
                return
            if not exists(path):
                savenpy(path, self.merge(parts))
            #if not exists(path):
            #    for key in inversion_set:
            #        unix.mkdir(path +'/'+ key)
            #        for proc in range(PAR.NPROC):
            #            with open(path +'/'+ key +'/'+ '%06d' % proc, 'w') as file:
            #                np.save(file, parts[key][proc])


    ### input file writers

    def write_parameters(self):
        unix.cd(self.getpath)
        write_parameters(vars(PAR))

    def write_receivers(self):
        unix.cd(self.getpath)
        key = 'use_existing_STATIONS'
        val = '.true.'
        solvertools.setpar(key, val)
        _, h = preprocess.load('traces/obs')
        write_receivers(h.nr, h.rx, h.rz)

    def write_sources(self):
        unix.cd(self.getpath)
        _, h = preprocess.load(dir='traces/obs')
        write_sources(vars(PAR), h)


    ### utility functions

    def mpirun(self, script, output='/dev/null'):
        """ Wrapper for mpirun
        """
        with open(output,'w') as f:
            subprocess.call(
                script,
                shell=True,
                stdout=f)

    @property
    def getlist(self):
        """list of all sources"""
        try:
            return self.events
        except:
            paths = glob(PATH.SOLVER_FILES +'/'+ 'SOURCE_*')
            names = [unix.basename(name) for name in paths]
            events = [name.split('_')[-1] for name in names]
            events.sort()
            self.events = events
            return self.events

    @ property
    def getname(self):
        """name of current source"""
        return self.getlist[system.getnode()]

    @property
    def getpath(self):
        """path of current source"""
        return join(PATH.SOLVER, self.getname)

    @property
    def wildcard(self):
        return glob('OUTPUT_FILES/U?_file_single.su')

    @property
    def databases(self):
        return join(self.getpath, 'OUTPUT_FILES/DATABASES_MPI')


