
import subprocess
from glob import glob
from os.path import join

import numpy as np

import seisflows.seistools.specfem3d as solvertools
from seisflows.seistools.shared import getpar, setpar

from seisflows.tools import unix
from seisflows.tools.code import exists, mpicall
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system


class specfem3d(custom_import('solver', 'base')):
    """ Python interface for SPECFEM3D

      See base class for method descriptions
    """

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
        mpicall(system.mpiexec(), 'bin/xspecfem3D')

        unix.mv(self.getdata, 'traces/obs')
        self.export_traces(PATH.OUTPUT, 'traces/obs')


    def generate_mesh(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(model_type)

        self.initialize_solver_directories()
        unix.cd(self.getpath)

        if model_type in ['gll']:
            par = getpar('MODEL').strip()
            if par != 'gll':
                if self.getnode == 0:
                    print 'WARNING: Unexpected Par_file setting:'
                    print 'MODEL =', par
            
            assert(exists(model_path))
            self.check_mesh_properties(model_path)

            src = glob(model_path +'/'+ '*')
            dst = self.model_databases
            unix.cp(src, dst)

            mpicall(system.mpiexec(), 'bin/xmeshfem3D')
            mpicall(system.mpiexec(), 'bin/xgenerate_databases')
            self.export_model(PATH.OUTPUT +'/'+ model_name)

        else:
            raise NotImplementedError


    ### low-level solver interface

    def forward(self):
        """ Calls SPECFEM3D forward solver
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        mpicall(system.mpiexec(), 'bin/xgenerate_databases')
        mpicall(system.mpiexec(), 'bin/xspecfem3D')


    def adjoint(self):
        """ Calls SPECFEM3D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')
        mpicall(system.mpiexec(), 'bin/xspecfem3D')

        # work around SPECFEM3D conflicting name conventions
        self.rename_data()


    ### input file writers

    def check_solver_parameter_files(self):
        """ Checks solver parameters
        """
        nt = getpar('NSTEP', cast=int)
        dt = getpar('DT', cast=float)

        if nt != PAR.NT:
            if self.getnode == 0: print "WARNING: nt != PAR.NT"
            setpar('NSTEP', PAR.NT)

        if dt != PAR.DT:
            if self.getnode == 0: print "WARNING: dt != PAR.DT"
            setpar('DT', PAR.DT)

        if self.mesh.nproc != PAR.NPROC:
            if self.getnode == 0:
                print 'Warning: mesh.nproc != PAR.NPROC'

        if 'MULTIPLES' in PAR:
            raise NotImplementedError


    def initialize_adjoint_traces(self):
        """ Works around SPECFEM3D file format issue by overriding base method
        """

    def initialize_adjoint_traces(self):
        super(specfem3d, self).initialize_adjoint_traces()
        self.rename_data()        

        # hack to deal with SPECFEM3D's requirement that all components exist,
        # even ones not in use
        unix.cd(self.getpath +'/'+ 'traces/adj')
        for iproc in range(PAR.NPROC):
            for channel in ['x', 'y', 'z']:
                src = '%d_d%s_SU.adj' % (iproc, PAR.CHANNELS[0])
                dst = '%d_d%s_SU.adj' % (iproc, channel)
                if not exists(dst):
                    unix.cp(src, dst)

    def rename_data(self):
        """ Works around conflicting data filename conventions
        """
        files = glob(self.getpath +'/'+ 'traces/adj/*SU')
        unix.rename('_SU', '_SU.adj', files)


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
        channels = PAR.CHANNELS
        return '*_d[%s]_SU' % channels.lower()

    @property
    def kernel_databases(self):
        return join(self.getpath, 'OUTPUT_FILES/DATABASES_MPI')

    @property
    def model_databases(self):
        return join(self.getpath, 'OUTPUT_FILES/DATABASES_MPI')

    @property
    def source_prefix(self):
        return 'FORCESOLUTION'

