
import subprocess
from os.path import join
from glob import glob

import numpy as np

import seisflows.seistools.specfem2d as solvertools
from seisflows.seistools.shared import getpar, setpar
from seisflows.seistools.io import splitvec, loadbypar

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import preprocess


class specfem2d(loadclass('solver', 'base')):
    """ Python interface for SPECFEM2D

      See base class for method descriptions
    """
    if PAR.MATERIALS == 'LegacyAcoustic':
        parameters = []
        parameters += ['vs']


    def check(self):
        """ Checks parameters and paths
        """
        super(specfem2d, self).check()

        # check time stepping parameters
        if 'NT' not in PAR:
            raise Exception

        if 'DT' not in PAR:
            raise Exception

        if 'F0' not in PAR:
            raise Exception


    def check_solver_parameter_files(self):
        """ Checks solver parameters
        """
        nt = getpar('nt', cast=int)
        dt = getpar('deltat', cast=float)
        f0 = getpar('f0', file='DATA/SOURCE', cast=float)

        if nt != PAR.NT:
            if system.getnode() == 0: print "WARNING: nt != PAR.NT"
            setpar('nt', PAR.NT)

        if dt != PAR.DT:
            if system.getnode() == 0: print "WARNING: dt != PAR.DT"
            setpar('deltat', PAR.DT)

        if f0 != PAR.F0:
            if system.getnode() == 0: print "WARNING: f0 != PAR.F0"
            setpar('f0', PAR.F0, file='DATA/SOURCE')

        if self.mesh.nproc != PAR.NPROC:
            if system.getnode() == 0:
                print 'Warning: mesh.nproc != PAR.NPROC'

        if 'MULTIPLES' in PAR:
            if PAR.MULTIPLES:
                setpar('absorbtop', '.false.')
            else:
                setpar('absorbtop', '.true.')


    def generate_data(self, **model_kwargs):
        """ Generates data
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.getpath)
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D',output='log.solver')

        unix.mv(self.data_wildcard, 'traces/obs')
        self.export_traces(PATH.OUTPUT, 'traces/obs')


    def generate_mesh(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(model_type)

        self.initialize_solver_directories()
        unix.cd(self.getpath)

        assert(exists(model_path))
        self.check_mesh_properties(model_path)

        src = glob(join(model_path, '*'))
        dst = join(self.getpath, 'DATA')
        unix.cp(src, dst)

        self.export_model(PATH.OUTPUT +'/'+ model_name)


    ### low-level solver interface

    def forward(self):
        """ Calls SPECFEM2D forward solver
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D')


    def adjoint(self):
        """ Calls SPECFEM2D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')

        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D')


    ### postprocessing utilities

    def smooth(self, path='', span=0., parameters=[]):
        """ Smooths SPECFEM2D kernels by convolving them with a Gaussian
        """
        from seisflows.tools.array import meshsmooth, stack

        assert parameters != []

        # implementing nproc > 1 would be straightforward, but a bit tedious
        assert self.mesh.nproc == 1

        kernels = self.load(path, suffix='_kernel')
        if not span:
            return kernels

        # set up grid
        _,x = loadbypar(PATH.MODEL_INIT, ['x'], 0)
        _,z = loadbypar(PATH.MODEL_INIT, ['z'], 0)
        mesh = stack(x[0], z[0])

        for key in self.parameters:
            kernels[key] = [meshsmooth(kernels[key][0], mesh, span)]

        unix.mv(path, path + '_nosmooth')
        self.save(path, kernels, suffix='_kernel')


    ### file transfer utilities

    def import_model(self, path):
        src = glob(path +'/'+ 'model/*')
        dst = join(self.getpath, 'DATA/')
        unix.cp(src, dst)

    def export_model(self, path):
        if system.getnode() == 0:
            unix.mkdir(path)
            src = glob(join(self.getpath, 'DATA/*.bin'))
            dst = path
            unix.cp(src, dst)


    ### input file writers

    def write_parameters(self):
        unix.cd(self.getpath)
        solvertools.write_parameters(vars(PAR))

    def write_receivers(self):
        unix.cd(self.getpath)
        key = 'use_existing_STATIONS'
        val = '.true.'
        setpar(key, val)
        _, h = preprocess.load('traces/obs')
        solvertools.write_receivers(h.nr, h.rx, h.rz)

    def write_sources(self):
        unix.cd(self.getpath)
        _, h = preprocess.load(dir='traces/obs')
        solvertools.write_sources(vars(PAR), h)


    ### utility functions

    def mpirun(self, script, output='/dev/null'):
        """ Wrapper for mpirun
        """
        with open(output,'w') as f:
            subprocess.call(
                script,
                shell=True,
                stdout=f)

    ### miscellaneous

    @property
    def data_wildcard(self):
        return glob('OUTPUT_FILES/U?_file_single.su')

    @property
    def model_databases(self):
        return join(self.getpath, 'DATA')

    @property
    def kernel_databases(self):
        return join(self.getpath, 'OUTPUT_FILES')

    @property
    def source_prefix(self):
        return 'SOURCE'


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
            self.mpirun(
                PATH.SPECFEM_BIN +'/'+ 'xcombine_sem '
                + name + '_kernel' + ' '
                + 'kernel_paths' + ' '
                + path +'/'+ 'sum')


    def generate_precond(self, process_traces=None, model_path=None, model_name=None, model_type='gll'):
        assert(model_name)
        assert(model_type)
        assert (exists(model_path))

        self.initialize_solver_directories()
        unix.cd(self.getpath)

        assert(exists(model_path))
        self.check_mesh_properties(model_path)

        src = glob(join(model_path, '*'))
        dst = join(self.getpath, 'DATA')
        unix.cp(src, dst)

        self.export_model(PATH.OUTPUT +'/'+ model_name)

        self.forward()
        unix.mv(self.data_wildcard, 'traces/syn')
        self.initialize_adjoint_traces('traces/syn')
        process_traces(self.getpath)

        self.adjoint()
        self.export_kernels(PATH.GLOBAL)

