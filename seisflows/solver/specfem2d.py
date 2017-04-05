
import sys
from os.path import basename, join
from glob import glob

import numpy as np

from seisflows.plugins.io import sem
from seisflows.tools.shared import getpar, setpar

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.tools import exists, call_solver, call_solver_nompi
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']


class specfem2d(custom_import('solver', 'base')):
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

        if 'WITH_MPI' not in PAR:
            setattr(PAR, 'WITH_MPI', False)

        # check time stepping parameters
        if 'NT' not in PAR:
            raise Exception

        if 'DT' not in PAR:
            raise Exception

        if 'F0' not in PAR:
            raise Exception

        # check data format
        if 'FORMAT' not in PAR:
            raise Exception()

        if PAR.FORMAT != 'su':
            raise Exception()


    def check_solver_parameter_files(self):
        """ Checks solver parameters
        """
        nt = getpar('nt', cast=int)
        dt = getpar('deltat', cast=float)
        f0 = getpar('f0', file='DATA/SOURCE', cast=float)

        if nt != PAR.NT:
            if self.getnode == 0: print "WARNING: nt != PAR.NT"
            setpar('nt', PAR.NT)

        if dt != PAR.DT:
            if self.getnode == 0: print "WARNING: dt != PAR.DT"
            setpar('deltat', PAR.DT)

        if f0 != PAR.F0:
            if self.getnode == 0: print "WARNING: f0 != PAR.F0"
            setpar('f0', PAR.F0, file='DATA/SOURCE')

        if self.mesh_properties.nproc != PAR.NPROC:
            if self.getnode == 0:
                print 'Warning: mesh_properties.nproc != PAR.NPROC'

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
        setpar('SAVE_FORWARD', '.false.')

        if PAR.WITH_MPI:
            call_solver(system.mpiexec(), 'bin/xmeshfem2D')
            call_solver(system.mpiexec(), 'bin/xspecfem2D')
        else:
            call_solver_nompi('bin/xmeshfem2D')
            call_solver_nompi('bin/xspecfem2D')

        if PAR.FORMAT in ['SU', 'su']:
            src = glob('OUTPUT_FILES/*.su')
            dst = 'traces/obs'
            unix.mv(src, dst)


    def initialize_adjoint_traces(self):
        super(specfem2d, self).initialize_adjoint_traces()

        # work around SPECFEM2D's use of different name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT in ['SU', 'su']:
            files = glob('traces/adj/*.su')
            unix.rename('.su', '.su.adj', files)

        # work around SPECFEM2D's requirement that all components exist,
        # even ones not in use
        if PAR.FORMAT in ['SU', 'su']:
            unix.cd(self.getpath +'/'+ 'traces/adj')
            for channel in ['x', 'y', 'z', 'p']:
                src = 'U%s_file_single.su.adj' % PAR.CHANNELS[0]
                dst = 'U%s_file_single.su.adj' % channel
                if not exists(dst):
                    unix.cp(src, dst)


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

    def forward(self, path='traces/syn'):
        """ Calls SPECFEM2D forward solver
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')

        if PAR.WITH_MPI:
            call_solver(system.mpiexec(), 'bin/xmeshfem2D')
            call_solver(system.mpiexec(), 'bin/xspecfem2D')
        else:
            call_solver_nompi('bin/xmeshfem2D')
            call_solver_nompi('bin/xspecfem2D')

        if PAR.FORMAT in ['SU', 'su']:
            filenames = glob('OUTPUT_FILES/*.su')
            unix.mv(filenames, path)


    def adjoint(self):
        """ Calls SPECFEM2D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')

        # hack to deal with SPECFEM2D's use of different name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT in ['SU', 'su']:
            files = glob('traces/adj/*.su')
            unix.rename('.su', '.su.adj', files)

        if PAR.WITH_MPI:
            call_solver(system.mpiexec(), 'bin/xmeshfem2D')
            call_solver(system.mpiexec(), 'bin/xspecfem2D')
        else:
            call_solver_nompi('bin/xmeshfem2D')
            call_solver_nompi('bin/xspecfem2D')


    ### postprocessing utilities

    def smooth(self, path='', parameters=None, span=0. ):
        """ For a long time SPECFEM2D lacked its own smoothing utility; this 
          method was intended only as a crude workaround
        """
        from seisflows.tools import array

        assert self.mesh_properties.nproc == 1,\
            msg.SmoothingError_SPECFEM2D

        kernels = self.load(path, suffix='_kernel')
        if not span:
            return kernels

        # set up grid
        x = sem.read(PATH.MODEL_INIT, 'x', 0)
        z = sem.read(PATH.MODEL_INIT, 'z', 0)
        mesh = array.stack(x, z)

        for key in parameters or self.parameters:
            kernels[key] = [array.meshsmooth(kernels[key][0], mesh, span)]

        unix.rm(path + '_nosmooth')
        unix.mv(path, path + '_nosmooth')
        self.save(path, kernels, suffix='_kernel')


    ### file transfer utilities

    def import_model(self, path):
        src = glob(path +'/'+ 'model/*')
        dst = join(self.getpath, 'DATA/')
        unix.cp(src, dst)

    def export_model(self, path):
        if self.getnode == 0:
            unix.mkdir(path)
            src = glob(join(self.getpath, 'DATA/*.bin'))
            dst = path
            unix.cp(src, dst)


    @property
    def data_filenames(self):
        if PAR.CHANNELS:
            if PAR.FORMAT in ['SU', 'su']:
               filenames = []
               for channel in PAR.CHANNELS:
                   filenames += ['U%s_file_single.su' % channel]
               return filenames

        else:
            unix.cd(self.getpath)
            unix.cd('traces/obs')

            if PAR.FORMAT in ['SU', 'su']:
                return glob('U?_file_single.su')

    @property
    def model_databases(self):
        return join(self.getpath, 'DATA')

    @property
    def kernel_databases(self):
        return join(self.getpath, 'OUTPUT_FILES')

    @property
    def source_prefix(self):
        return 'SOURCE'

