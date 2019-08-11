#
# This is Seisflows
#
# See LICENCE file
#
###############################################################################

# Import system modules
import sys
from os.path import basename, join
from glob import glob

# Import Numpy
import numpy as np

# Local imports
from seisflows.plugins.solver.specfem2d import smooth_legacy
from seisflows.tools.seismic import getpar, setpar
from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.seismic import call_solver
from seisflows.tools.tools import exists
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
        found = False
        for ntStr in ['NSTEP', 'nt']:  # Different specfem conventions
            try:
                nt = getpar(ntStr, cast=int, noOutput=True)
                found = True
            except Exception:
                pass
        if not found:
            raise Exception('Not found in Specfem Par_file: NSTEP or nt')
        found = False
        for dtStr in ['DT', 'deltat']:  # Different specfem conventions
            try:
                dt = getpar(dtStr, cast=float, noOutput=True)
                found = True
            except Exception:
                pass
        if not found:
            raise Exception('Not found in Specfem Par_file: DT or deltat')

        f0 = getpar('f0', file='DATA/SOURCE', cast=float)

        if nt != PAR.NT:
            if self.taskid == 0:
                print "WARNING: nt != PAR.NT", nt, PAR.NT
            setpar('nt', PAR.NT)

        if dt != PAR.DT:
            if self.taskid == 0:
                print "WARNING: dt != PAR.DT", dt, PAR.DT
            setpar('deltat', PAR.DT)

        if f0 != PAR.F0:
            if self.taskid == 0:
                print "WARNING: f0 != PAR.F0"
            setpar('f0', PAR.F0, filename='DATA/SOURCE')

        if self.mesh_properties.nproc != PAR.NPROC:
            if self.taskid == 0:
                print 'Warning: mesh_properties.nproc != PAR.NPROC'

        if 'MULTIPLES' in PAR:
            if PAR.MULTIPLES:
                setpar('absorbtop', '.false.')
            else:
                setpar('absorbtop', '.true.')

    def generate_data(self, **model_kwargs):
        """ Generates data (perform meshing and database generation first)
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.cwd)
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.false.')

        call_solver(system.mpiexec(), 'bin/xmeshfem2D', output='mesher.log')
        call_solver(system.mpiexec(), 'bin/xspecfem2D', output='solver.log')

        if PAR.FORMAT in ['SU', 'su']:
            src = glob('OUTPUT_FILES/*.su')
            # work around SPECFEM2D's different file names (depending on the
            # version used :
            unix.rename('single_p.su', 'single.su', src)
            src = glob('OUTPUT_FILES/*.su')
            dst = 'traces/obs'
            unix.mv(src, dst)

        if PAR.SAVETRACES:
            self.export_traces(PATH.OUTPUT+'/'+'traces/obs')

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
            unix.cd(self.cwd + '/' + 'traces/adj')
            for channel in ['x', 'y', 'z', 'p']:
                src = 'U%s_file_single.su.adj' % PAR.CHANNELS[0]
                dst = 'U%s_file_single.su.adj' % channel
                if not exists(dst):
                    unix.cp(src, dst)

    def generate_mesh(self, model_path=None, model_name=None,
                      model_type='gll'):
        """ Performs meshing and database generation
        """
        # Assert that the model name and path are not empty
        assert(model_name)
        assert(model_type)

        self.initialize_solver_directories()
        unix.cd(self.cwd)

        assert(exists(model_path))  # Check that the path exist
        # Fill _mesh_properties which contain the number of integration points,
        # the number of procs used and the coordinates of the points
        self.check_mesh_properties(model_path)

        # Copy the model files (ex: proc000023_vp.bin ...) into DATA
        src = glob(join(model_path, '*'))
        dst = join(self.cwd, 'DATA')
        unix.cp(src, dst)

        # Export the model into output folder
        if self.taskid == 0:
            self.export_model(PATH.OUTPUT + '/' + model_name)

    # Low-level solver interface

    def forward(self, path='traces/syn'):
        """ Calls SPECFEM2D forward solver
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')

        call_solver(system.mpiexec(), 'bin/xmeshfem2D')
        call_solver(system.mpiexec(), 'bin/xspecfem2D')

        if PAR.FORMAT in ['SU', 'su']:
            filenames = glob('OUTPUT_FILES/*.su')
            # work around SPECFEM2D's different file names (depending on the
            # version used :
            unix.rename('single_p.su', 'single.su', filenames)
            filenames = glob('OUTPUT_FILES/*.su')
            unix.mv(filenames, path)

    def adjoint(self):
        """ Calls SPECFEM2D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')

        # hack to deal with different SPECFEM2D name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT in ['SU', 'su']:
            files = glob('traces/adj/*.su')
            unix.rename('.su', '.su.adj', files)

        call_solver(system.mpiexec(), 'bin/xmeshfem2D')
        call_solver(system.mpiexec(), 'bin/xspecfem2D')

    # File transfer utilities

    def import_model(self, path):
        src = glob(path + '/' + 'model/*')
        dst = join(self.cwd, 'DATA/')
        unix.cp(src, dst)

    def export_model(self, path):
        unix.mkdir(path)
        src = glob(join(self.cwd, 'DATA/*.bin'))
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
            unix.cd(self.cwd)
            unix.cd('traces/obs')

            if PAR.FORMAT in ['SU', 'su']:
                return glob('U?_file_single.su')

    @property
    def model_databases(self):
        return join(self.cwd, 'DATA')

    @property
    def kernel_databases(self):
        return join(self.cwd, 'OUTPUT_FILES')

    @property
    def source_prefix(self):
        return 'SOURCE'

    # workaround for older versions of SPECFEM2D,
    # which lacked a smoothing utility
    # if not exists(PATH.SPECFEM_BIN+'/'+'xsmooth_sem'):
    #    smooth = staticmethod(smooth_legacy)
    smooth = staticmethod(smooth_legacy)
