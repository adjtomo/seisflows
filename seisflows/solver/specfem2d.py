
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
        assert (exists(model_path))

        self.initialize_solver_directories()

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

    def smooth(self, path='', span=0.):
        """ Smooths SPECFEM2D kernels by convolving them with a Gaussian
        """
        from seisflows.tools.array import meshsmooth, stack

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
        return join(self.getpath, 'OUTPUT_FILES')

    @property
    def source_prefix(self):
        return 'SOURCE'

