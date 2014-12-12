import subprocess

import numpy as np

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.code import exists, glob, join
from seisflows.tools.config import findpath, ConfigObj, ParameterObj
from seisflows.tools.io import loadbin, savebin

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class specfem3d_globe(object):
    """ Python interface for SPECFEM3D_GLOBE

      eval_func, eval_grad, apply_hess
        These methods deal with evaluation of the misfit function or its
        derivatives and provide the primary interface between the solver and
        other workflow components.

      forward, adjoint, mesher
        These methods allow direct access to individual SPECFEM3D_GLOBE
        components. Together, they provide a secondary interface users can
        employ for specialized tasks not covered by high level methods.

      prepare_solver, prepare_data, prepare_model
        SPECFEM3D_GLOBE requires a particular directory structure in which to
        run and particular file formats for models, data, and parameter files.
        These methods help put in place all these prerequisites.

      load, save
        For reading and writing SPECFEM3D_GLOBE models and kernels. On the disk,
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

    if 0:
        # isotropic
        model_parameters = []
        model_parameters += ['reg1_rho']
        model_parameters += ['reg1_vp']
        model_parameters += ['reg1_vs']
        # model_parameters += ['reg2_rho']
        # model_parameters += ['reg2_vp']
        # model_parameters += ['reg2_vs']
        # model_parameters += ['reg3_rho']
        # model_parameters += ['reg3_vp']
        # model_parameters += ['reg3_vs']

        inversion_parameters = []
        inversion_parameters += ['reg1_rho']
        inversion_parameters += ['reg1_vp']
        inversion_parameters += ['reg1_vs']

        kernel_map = {
            'reg1_rho': 'reg1_rho_kernel',
            'reg1_vp': 'reg1_alpha_kernel',
            'reg1_vs': 'reg1_beta_kernel'}
        # 'reg2_rho':'reg2_rho_kernel',
        # 'reg2_vp':'reg2_alpha_kernel',
        # 'reg2_vs':'reg2_beta_kernel',
        # 'reg3_rho':'reg3_rho_kernel',
        # 'reg3_vp':'reg3_alpha_kernel',
        # 'reg3_vs':'reg3_beta_kernel'}

    else:
        # transversely isotropic
        model_parameters = []
        model_parameters += ['reg1_rho']
        model_parameters += ['reg1_vpv']
        model_parameters += ['reg1_vph']
        model_parameters += ['reg1_vsv']
        model_parameters += ['reg1_vsh']
        model_parameters += ['reg1_eta']
        # model_parameters += ['reg2_rho']
        # model_parameters += ['reg2_vp']
        # model_parameters += ['reg2_vs']
        # model_parameters += ['reg3_rho']
        # model_parameters += ['reg3_vp']
        # model_parameters += ['reg3_vs']

        inversion_parameters = []
        inversion_parameters += ['reg1_rho']
        inversion_parameters += ['reg1_vpv']
        inversion_parameters += ['reg1_vph']
        inversion_parameters += ['reg1_vsv']
        inversion_parameters += ['reg1_vsh']
        inversion_parameters += ['reg1_eta']

        kernel_map = {
            'reg1_rho': 'reg1_rho_kernel',
            'reg1_eta': 'reg1_rho_kernel',
            'reg1_vph': 'reg1_alpha_kernel',
            'reg1_vpv': 'reg1_alpha_kernel',
            'reg1_vsv': 'reg1_beta_kernel',
            'reg1_vsh': 'reg1_beta_kernel'}
        # 'reg2_rho':'reg2_rho_kernel',
        # 'reg2_vp':'reg2_alpha_kernel',
        # 'reg2_vs':'reg2_beta_kernel',
        # 'reg3_rho':'reg3_rho_kernel',
        # 'reg3_vp':'reg3_alpha_kernel',
        # 'reg3_vs':'reg3_beta_kernel'}

    # data channels
    channels = []
    channels += ['z']

    # data input/output
    reader = staticmethod(seistools.specfem3d_globe.read)
    writer = staticmethod(seistools.specfem3d_globe.write)
    glob = lambda _: glob('OUTPUT_FILES/*.sem.ascii')

    def check(self):
        """ Checks parameters, paths, and dependencies
        """

        # check parameters
        if 'NEX_XI' not in PAR:
            pass

        if 'NEX_ETA' not in PAR:
            pass

        if 'NPROC_XI' not in PAR:
            pass

        if 'NPROC_ETA' not in PAR:
            pass

        # check paths
        if 'GLOBAL' not in PATH:
            raise Exception

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'MESH' not in PATH:
            setattr(PATH, 'MESH', join(PATH.GLOBAL, 'mesh'))

        if 'SOLVER' not in PATH:
            if PATH.LOCAL:
                setattr(PATH, 'SOLVER', join(PATH.LOCAL, 'solver'))
            else:
                setattr(PATH, 'SOLVER', join(PATH.GLOBAL, 'solver'))

        # check dependencies
        if 'preprocess' not in OBJ:
            raise Exception

        if 'system' not in OBJ:
            raise Exception("Undefined Exception")

        global preprocess
        import preprocess

        global system
        import system

    def setup(self):
        """ Prepares solver for inversion, migration, or forward modeling
        """
        self.prepare_dirs()
        model_type = seistools.specfem3d_globe.getpar('MODEL')

        # prepare data
        if exists(PATH.DATA):
            self.prepare_data(
                data_path=PATH.DATA)
        else:
            self.prepare_data(
                model_path=PATH.MODEL_TRUE,
                model_name='model_true',
                model_type=model_type)

        # prepare model
        self.prepare_model(
            model_path=PATH.MODEL_INIT,
            model_name='model_init',
            model_type=model_type)

    def prepare_dirs(self):
        """ Sets up directory in which to run solver

          Creates subdirectories expected by SPECFEM3D, copies mesher and solver
          binary files, and optionally calls prepare_data and prepare_mdoel.
          Binaries must be supplied by user as there is currently no mechanism
          to automatically compile from source code.
        """
        unix.rm(self.path)
        unix.mkdir(self.path)
        unix.cd(self.path)

        # create subdirectories
        unix.mkdir('bin')
        unix.mkdir('DATA')
        unix.mkdir('OUTPUT_FILES/DATABASES_MPI')

        unix.mkdir('traces/obs')
        unix.mkdir('traces/syn')
        unix.mkdir('traces/adj')

        # copy binaries
        src = glob(PATH.SOLVER_BINARIES + '/' + '*')
        dst = 'bin/'
        unix.cp(src, dst)

        # copy input files
        src = glob(PATH.SOLVER_FILES + '/' + '*')
        dst = 'DATA/'
        unix.cp(src, dst)

    def prepare_data(self, data_path=None, **kwargs):
        """ Prepares data for inversion or migration

          Users implementing an inversion or migration can choose between
          supplying data or supplying a target model from which data are
          generated on the fly. In both cases, all necessary SPECFEM3D input
          files must be provided.

          Adjoint traces are intialized by writing zeros for all components.
          This is necessary because, at the start of an adjoint simulation,
          SPECFEM3D expects that all components exists, even ones not actually
          in use for the inversion.
        """
        unix.cd(self.path)

        # update source coordinates
        src = 'DATA/FORCESOLUTION_' + self.getshot()
        dst = 'DATA/FORCESOLUTION'
        unix.cp(src, dst)

        if data_path:
            # copy user supplied data
            src = glob(data_path + '/' + self.getshot() + '/' + '*')
            dst = 'traces/obs/'
            unix.cp(src, dst)
            self.initialize_adjoint()

        else:
            # generate data
            self.prepare_model(**kwargs)
            self.mpirun('bin/xspecfem3D')
            unix.mv(self.glob(), 'traces/obs')
            self.export_traces(PATH.OUTPUT, 'traces/obs')
            self.initialize_adjoint()

    def prepare_model(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and model interpolation using SPECFEM3D's builtin
          mesher and database generation utility
        """
        assert model_type
        print 'model_path:', model_path
        print 'model_name:', model_name
        print 'model_type:', model_type
        print ''

        unix.cd(self.path)

        # run builtin mesher and generate databases
        if model_type == 'gll':
            assert (exists(model_path))
            # copy files
            src = glob(model_path + '/' + '*')
            dst = 'OUTPUT_FILES/DATABASES_MPI/'
            unix.cp(src, dst)

        else:
            pass

        seistools.specfem3d_globe.setpar('MODEL', model_type)
        self.mesher()

        # save results
        parts = self.load('OUTPUT_FILES/DATABASES_MPI')
        if system.getnode() == 0:
            if model_name and model_type == 'gll':
                unix.ln(model_path, PATH.OUTPUT + '/' + model_name)
            elif model_name:
                self.save(PATH.OUTPUT + '/' + model_name, parts)
            if not exists(PATH.MESH):
                set1 = set(self.model_parameters)
                set2 = set(self.inversion_parameters)
                keys = list(set1.difference(set2))
                for key in keys:
                    unix.mkdir(PATH.MESH + '/' + key)
                    for proc in range(PAR.NPROC):
                        with open(PATH.MESH + '/' + key + '/' + '%06d' % proc,
                                  'w') as file:
                            np.save(file, parts[key][proc])

    # -- high-level solver interface

    def eval_func(self, path='', export_traces=False):
        """ Evaluates misfit function by carrying out forward simulation and
            making measurements on observations and synthetics.
        """
        unix.cd(self.path)
        self.import_model(path)

        # forward simulation
        self.forward()
        unix.mv(self.glob(), 'traces/syn')
        preprocess.prepare_eval_grad(self.path)

        # save results
        self.export_residuals(path)
        if export_traces:
            self.export_traces(path, prefix='traces/syn')

    def eval_grad(self, path='', export_traces=False):
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
            self.export_traces(path, prefix='traces/syn')

    def apply_hess(self, path='', hessian='Newton'):
        """ Evaluates action of Hessian on a given model vector.
        """
        unix.cd(self.path)
        self.imprt(path, 'model')

        # forward simulation
        self.forward()
        unix.mv(self.glob(), 'traces/lcg')

        preprocess.prepare_apply_hess(self.path)

        # adjoint simulation
        self.adjoint()

        # save results
        # FIXME: declaration has only 1 argument. Should it have 2?
        self.export_kernels(path, 'kernels')

    # -- low-level solver interface

    def forward(self):
        """ Calls SPECFEM3D forward solver
        """
        # prepare solver
        seistools.specfem3d_globe.setpar('SIMULATION_TYPE', '1')
        seistools.specfem3d_globe.setpar('SAVE_FORWARD', '.true.')
        seistools.specfem3d_globe.setpar('MODEL', 'gll')

        # run solver
        self.mpirun('bin/xspecfem3D')
        unix.mv(self.glob(), 'traces/syn')

    def adjoint(self):
        """ Calls SPECFEM3D adjoint solver
        """
        # prepare solver
        seistools.specfem3d_globe.setpar('SIMULATION_TYPE', '3')
        seistools.specfem3d_globe.setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')

        # run solver
        self.mpirun('bin/xspecfem3D')

    def mesher(self):
        """ Calls SPECFEM3D builtin mesher
        """
        self.mpirun('bin/xmeshfem3D')

    # -- model input/output

    def load(self, dirname, type='model'):
        """ reads SPECFEM3D kernel or model to dictionary
        """
        mapping = lambda key: self.kernel_map[key]
        parts = {}

        if type == 'model':
            for key in self.model_parameters:
                parts[key] = []
                # read database files
                for iproc in range(PAR.NPROC):
                    filename = 'proc%06d_%s.bin' % (iproc, key)
                    part = loadbin(join(dirname, filename))
                    parts[key].append(part)

        elif type == 'kernel':
            for key in self.model_parameters:
                parts[key] = []
                # read database files
                for iproc in range(PAR.NPROC):
                    filename = 'proc%06d_%s.bin' % (iproc, mapping(key))
                    part = loadbin(join(dirname, filename))
                    parts[key].append(part)

        return parts

    def save(self, dirname, parts):
        """ writes SPECFEM3D model
        """
        unix.mkdir(dirname)

        # write database files
        for key in self.model_parameters:
            nn = len(parts[key])
            for ii in range(nn):
                filename = 'proc%06d_%s.bin' % (ii, key)
                savebin(parts[key][ii], join(dirname, filename))

    # -- vector/dictionary conversion

    def merge(self, parts):
        """ merges dictionary into vector
        """
        v = np.array([])
        for key in self.inversion_parameters:
            for iproc in range(PAR.NPROC):
                v = np.append(v, parts[key][iproc])

    def split(self, v):
        """ splits vector into dictionary
        """
        parts = {}
        nrow = len(v)/(PAR.NPROC*len(self.inversion_parameters))
        j = 0
        for key in self.model_parameters:
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
                    parts[key].append(
                        np.load(PATH.MESH + '/' + key + '/' + proc))
        return parts

    # -- postprocessing utilities

    def combine(self, path=''):
        """ combines SPECFEM3D kernels
        """
        dirs = unix.ls(path)

        # initialize kernels
        unix.mkdir(path + '/' + 'sum')
        for key in self.model_parameters:
            if key not in self.inversion_parameters:
                for i in range(PAR.NPROC):
                    proc = '%06d' % i
                    src = PATH.MESH + '/' + key + '/' + proc
                    dst = path + '/' + 'sum' + '/' + 'proc' + proc + '_' + \
                        self.kernel_map[key] + '.bin'
                    savebin(np.load(src), dst)

        # create temporary files and directories needed by xsum_kernels
        unix.cd(self.path)
        with open('kernels_list.txt', 'w') as file:
            file.write('\n'.join(dirs) + '\n')
        unix.mkdir('INPUT_KERNELS')
        unix.mkdir('OUTPUT_SUM')
        for dir in dirs:
            src = path + '/' + dir
            dst = unix.pwd() + '/' + 'INPUT_KERNELS' + '/' + dir
            unix.ln(src, dst)

        # sum kernels
        self.mpirun(PATH.SOLVER_BINARIES + '/' + 'xsum_kernels')
        unix.mv(glob('OUTPUT_SUM/*'), path + '/' + 'sum')

        # remove temporary files and directories
        unix.rm('INPUT_KERNELS')
        unix.rm('OUTPUT_SUM')
        unix.rm('kernels_list.txt')
        unix.cd(path)

    def smooth(self, path='', tag='grad', span=0):
        """ smooths SPECFEM3D kernels
        """
        unix.mv(path + '/' + tag, path + '/' + tag + '_nosmooth')
        unix.mkdir(path + '/' + tag)

        # prepare list
        kernel_list = []
        for key in self.model_parameters:
            if key in self.inversion_parameters:
                smoothing_flag = True
            else:
                smoothing_flag = False
            region, kernel_name = key.split('_')
            kernel_list = kernel_list + [[kernel_name, smoothing_flag]]

        unix.cd(self.path)

        for kernel_name, smoothing_flag in kernel_list:
            if smoothing_flag:
                # run smoothing
                print ' smoothing', kernel_name
                self.mpirun(
                    PATH.SOLVER_BINARIES + '/' + 'xsmooth_sem '
                    + str(span) + ' '
                    + str(span) + ' '
                    + kernel_name + ' '
                    + path + '/' + tag + '_nosmooth/' + ' '
                    + self.path + '/'
                    + 'OUTPUT_FILES/DATABASES_MPI' + '/' + ' ')

                src = glob(path + '/' + tag + '_nosmooth' + '/*_smooth.bin')
                dst = path + '/' + tag
                unix.mv(src, dst)
                unix.rename('_smooth', '',
                            glob(path + '/' + tag + '/*_smooth.bin'))

            else:
                src = glob(
                    path + '/' + tag + '_nosmooth/*' + kernel_name + '.bin')
                dst = path + '/' + tag + '/'
                unix.cp(src, dst)

        print ''

        unix.cd(path)

    # -- input file writers

    def write_parameters(self):
        unix.cd(self.path)

        write_parameters(vars(PAR))

    def write_receivers(self):
        unix.cd(self.path)

        # adjust parameters
        key = 'use_existing_STATIONS'
        val = '.true.'
        seistools.specfem3d_globe.setpar(key, val)

        # write receivers file
        _, h = preprocess.load('traces/obs')
        write_receivers(h.nr, h.rx, h.rz)

    def write_sources(self):
        unix.cd(self.path)

        # write source file
        _, h = preprocess.load(dir='traces/obs')
        write_sources(vars(PAR), h)

    # -- file transfer utilities

    def import_model(self, path):
        src = glob(join(path, 'model', '*'))
        dst = join(unix.pwd(), 'OUTPUT_FILES/DATABASES_MPI')
        unix.cp(src, dst)

    def import_traces(self, path):
        src = glob(join(path, 'traces', getname(), '*'))
        dst = join(unix.pwd(), 'traces/obs')
        unix.cp(src, dst)

    def export_model(self, path):
        if system.getnode() != 0:
            return
        for name in self.model_parameters:
            src = glob(join(unix.pwd(), 'OUTPUT_FILES/DATABASES_MPI',
                            '*_' + name + '.bin'))
            dst = path
            unix.mkdir(dst)
            unix.cp(src, dst)

    def export_kernels(self, path):
        try:
            unix.mkdir(join(path, 'kernels'))
        except OSError:
            pass
        unix.mkdir(join(path, 'kernels', '%06d' % system.getnode()))
        for name in self.kernel_map.values():
            src = join(glob(
                unix.pwd() + '/' + 'OUTPUT_FILES/DATABASES_MPI'
                + '/' + '*' + name + '.bin'))
            dst = join(path, 'kernels', '%06d' % system.getnode())
            unix.mv(src, dst)
        try:
            name = 'rhop_kernel'
            src = join(glob(
                unix.pwd() + '/' + 'OUTPUT_FILES/DATABASES_MPI'
                + '/' + '*' + name + '.bin'))
            dst = join(path, 'kernels', '%06d' % system.getnode())
            unix.mv(src, dst)
        except OSError:
            pass

    def export_residuals(self, path):
        try:
            unix.mkdir(join(path, 'residuals'))
        except OSError:
            pass
        src = join(unix.pwd(), 'residuals')
        dst = join(path, 'residuals', '%06d' % system.getnode())
        unix.mv(src, dst)

    def export_traces(self, path, prefix='traces/obs'):
        try:
            unix.mkdir(join(path, 'traces'))
        except OSError:
            pass
        src = join(unix.pwd(), prefix)
        dst = join(path, 'traces', '%06d' % system.getnode())
        unix.cp(src, dst)

    def initialize_adjoint(self):
        _, h = preprocess.load('traces/obs')
        zeros = np.zeros((h.nt, h.nr))
        for channel in ['x', 'y', 'z']:
            self.writer(zeros, h, channel=channel, prefix='traces/adj')

    # -- utility functions

    def cleanup(self):
        """ Cleans up directory after simulation
        """
        unix.cd(self.path)
        unix.rm(glob('traces/syn/*'))
        unix.rm(glob('traces/adj/*'))

    def mpirun(self, script, outfile='/dev/null'):
        """ Wrapper for mpirun
        """
        with open(outfile) as f:
            subprocess.call(
                system.mpiargs()
                + script,
                shell=True,
                stdout=f)

    @property
    def path(self):
        return join(PATH.SOLVER, self.getshot())

    def getshot(self):
        return '%06d' % system.getnode()

    def gettype(self):
        try:
            return seistools.specfem3d_globe.getpar(
                'MODEL', file=PATH.SOLVER_FILES + '/' + 'Par_file')
        except:
            return None

