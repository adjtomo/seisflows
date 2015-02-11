
import subprocess
from glob import glob
from os.path import join

import numpy as np

import seisflows.seistools.specfem3d as solvertools
from seisflows.seistools.io import load
from seisflows.seistools.shared import getpar, setpar

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class specfem3d(loadclass('solver', 'base')):
    """ Python interface for SPECFEM3D

      See base class for method descriptions
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

    kernel_map = {
        'rho': 'rho_kernel',
        'vp': 'alpha_kernel',
        'vs': 'beta_kernel'}


    def check(self):
        """ Checks parameters and paths
        """
        super(specfem3d, self).check()

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
        self.mpirun('bin/xspecfem3D')

        unix.mv(self.data_wildcard, 'traces/obs')
        self.export_traces(PATH.OUTPUT, 'traces/obs')


    def generate_mesh(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(model_type)

        self.initialize_solver_directories()
        unix.cd(self.getpath)

        if model_type in ['gll']:
            assert (exists(model_path))
            unix.cp(glob(model_path +'/'+ '*'), self.model_databases)
            self.mpirun('bin/xmeshfem3D')
            self.mpirun('bin/xgenerate_databases')
            self.export_model(PATH.OUTPUT +'/'+ model_name)

        elif model_type in ['sep']:
            self.mpirun('bin/xmeshfem3D')
            self.mpirun('bin/xgenerate_databases')
            self.export_model(PATH.OUTPUT +'/'+ model_name)

        elif model_type == 'tomo':
            raise NotImplementedError


    ### low-level solver interface

    def forward(self):
        """ Calls SPECFEM3D forward solver
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xgenerate_databases')
        self.mpirun('bin/xspecfem3D')


    def adjoint(self):
        """ Calls SPECFEM3D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')
        self.mpirun('bin/xspecfem3D')


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


    ### miscellaneous

    @property
    def data_wildcard(self):
        return glob('OUTPUT_FILES/*SU')

    @property
    def model_databases(self):
        return join(self.getpath, 'OUTPUT_FILES/DATABASES_MPI')

    @property
    def source_prefix(self):
        return 'FORCESOLUTION'



