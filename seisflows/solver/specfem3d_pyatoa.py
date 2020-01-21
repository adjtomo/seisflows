#!/usr/bin/env python
"""
This is the subclass seisflows.solver.Specfem3D_pyatoa
This class provides utilities for the Seisflows solver interactions with
Specfem3D Cartesian. It inherits all attributes from seisflows.solver.Base,
and overwrites these functions to provide specified interaction with Specfem3D

This subclass also interacts with the Python package Pyatoa, to perform data
collection, preprocessing, and misfit quantification steps
"""
import os
import sys
import time
import subprocess
import seisflows.plugins.solver.specfem3d as solvertools

from glob import glob
from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import custom_import
from seisflows.tools.seismic import call_solver, getpar, setpar

from pyatoa.plugins.seisflows.pyaflowa import Pyaflowa

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']


class Specfem3dPyatoa(custom_import('solver', 'Base')):
    """
    Python interface to Specfem3D Cartesian. This subclass inherits functions
    from seisflows.solver.Base

    !!! See base class for method descriptions !!!
    """
    def check(self):
        """
        Checks parameters and paths. Most prechecks will happen in the workflow.
        """
        # Run Base class checks
        super(Specfem3dPyatoa, self).check()

        # Check time stepping parameters
        if "NT" not in PAR:
            raise Exception("'NT' not specified in parameters file")

        if "DT" not in PAR:
            raise Exception("'DT' not specified in parameters file")

        # Check data format for Specfem3D
        if "FORMAT" not in PAR:
            raise Exception("'FORMAT' not specified in parameters file")

        # Currently Pyatoa only works with ASCII
        if PAR.FORMAT not in ["ascii"]:
            raise Exception()

    def setup(self):
        """
        Overload of solver.base.setup

        Removes the need to move data around as Pyatoa takes care of data
        fetching within eval_func.

        Allows for synthetic-synthetic cases
        """
        # Clean up for new inversion
        unix.rm(self.cwd)
        self.initialize_solver_directories()
        
        # If synthetic case, create synthetic observations
        if PAR.CASE == "Synthetic":
            self.generate_data(model_path=PATH.MODEL_TRUE,
                               model_name="model_true", model_type="gll")

        # Prepare initial model
        self.generate_mesh(model_path=PATH.MODEL_INIT, model_name="model_init",
                           model_type="gll")

    def generate_data(self, **model_kwargs):
        """
        Overload seisflows.solver.base.generate_data
        
        Not used if PAR.CASE == "Data"

        Generates data in the synthetic-synthetic comparison case.
        Automatically calls generate mesh for the true model, rather than
        passing them in as kwargs.

        Also turns on attenuation for the forward model
        !!! attenuation could be moved into parameters.yaml? !!!
        """
        # Prepare for the forward simulation
        self.generate_mesh(**model_kwargs)
        if PAR.VERBOSE:
            print("Specfem3dPyatoa.generate data")

        unix.cd(self.cwd)
        setpar("SIMULATION_TYPE", "1")
        setpar("SAVE_FORWARD", ".true.")
        setpar("ATTENUATION ", ".true.")
        call_solver(mpiexec=system.mpiexec(), exectuable="bin/xspecfem3D")

        # ASCII .sem? output format
        if PAR.FORMAT == "ascii":
            unix.mv(src=glob(os.path.join("OUTPUT_FILES", "*sem?")),
                    dst="traces/obs")
        else:
            raise NotImplementedError("Pyatoa only works with ascii")

        # Export traces to permanent storage on disk
        if PAR.SAVETRACES:
            self.export_traces(os.path.join(PATH.OUTPUT, "traces", "obs"))

    def generate_mesh(self, model_path, model_name, model_type='gll'):
        """
        Performs meshing and database generation

        :type model_path: str
        :param model_path: path to the model to be used for mesh generation
        :type model_name: str
        :param model_name: name of the model to be used as identification
        :type model_type: str
        :param model_type: available model types to be passed to the Specfem3D
            Par_file. See Specfem3D Par_file for available options.
        """
        if PAR.VERBOSE:
            print("Specfem3dPyatoa.generate mesh")

        unix.cd(self.cwd)

        # Run mesh generation
        assert(exists(model_path))
        if model_type in available_model_types:
            par = getpar("MODEL").strip()
            if par == "gll":
                self.check_mesh_properties(model_path)

                src = glob(os.path.join(model_path, "*"))
                dst = self.model_databases
                unix.cp(src, dst)

                call_solver(mpiexec=system.mpiexec(),
                            executable="bin/xgenerate_databases")

            # Export the model for future use in the workflow
            if self.taskid == 0:
                self.export_model(os.path.join(PATH.OUTPUT, model_name))
        else:
            raise NotImplementedError(f"MODEL={par} not implemented")

    def eval_fwd(self, path=''):
        """
        Performs forward simulations for misfit function evaluation.
        Same as solver.base.eval_func without the residual writing.

        Function evaluation is taken care of by Pyatoa.
        """
        if PAR.VERBOSE:
            print("Specfem3dPyatoa.eval_fwd")
        unix.cd(self.cwd)
        self.import_model(path)
        self.forward()
    
    def eval_func(self, iter='', step=0, suffix=None, *args, **kwargs):
        """
        Call Pyatoa workflow to evaluate the misfit functional.

        :type iter:
        :param args:
        :param kwargs:
        :return:
        """
        if PAR.VERBOSE:
            print("Specfem3dPyatoa.eval_func calling Pyatoa")


        arguments = " ".join([
            "--mode process",
            "--working_dir {}".format(PATH.WORKDIR),
            "--current_dir {}".format(self.cwd),
            "--model_number {}".format("m{:0>2}".format(int(iter) - 1)),
            "--event_id {}".format(self.source_name),
            "--step_count {}".format("s{:0>2}".format(step)),
            "--suffix {}".format(suffix)
        ])
        call_pyatoa = " ".join([load_conda, load_hdf5, PATH.PYTHON3,
                                PATH.PYATOA_RUN, arguments
                                ])
        print(call_pyatoa)
        try:
            tstart = time.time()
            stdout = subprocess.check_output(call_pyatoa, shell=True)
            print '{:.2f}m elapsed'.format((time.time() - tstart) / 60 )
        except subprocess.CalledProcessError as e:
            print("Pyatoa failed with {}".format(e))
            sys.exit(-1)
    
    # low-level solver interface
    def forward(self, path='traces/syn'):
        """ Calls SPECFEM3D forward solver and then moves files into path
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        setpar('ATTENUATION ', '.true.')
        call_solver(system.mpiexec(), 'bin/xgenerate_databases')
        call_solver(system.mpiexec(), 'bin/xspecfem3D')

        # seismic unix output format
        if PAR.FORMAT in ['SU', 'su']:
            src = glob('OUTPUT_FILES/*_d?_SU')
            dst = path
            unix.mv(src, dst)
        # ascii sem output format
        elif PAR.FORMAT == "ascii":
            src = glob('OUTPUT_FILES/*sem?')
            dst = path
            unix.mv(src, dst)

    def adjoint(self):
        """ Calls SPECFEM3D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        setpar('ATTENUATION ', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')
        call_solver(system.mpiexec(), 'bin/xspecfem3D')

    # input file writers
    def check_solver_parameter_files(self):
        """ Checks solver parameters
        """
        nt = getpar('NSTEP', cast=int)
        dt = getpar('DT', cast=float)

        if nt != PAR.NT:
            if self.taskid == 0: print "WARNING: nt != PAR.NT"
            setpar('NSTEP', PAR.NT)

        if dt != PAR.DT:
            if self.taskid == 0: print "WARNING: dt != PAR.DT"
            setpar('DT', PAR.DT)

        if self.mesh_properties.nproc != PAR.NPROC:
            if self.taskid == 0:
                print 'Warning: mesh_properties.nproc != PAR.NPROC'

        if 'MULTIPLES' in PAR:
            raise NotImplementedError

    def initialize_adjoint_traces(self):
        super(Specfem3dPyatoa, self).initialize_adjoint_traces()

        # workaround for SPECFEM2D's use of different name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT in ['SU', 'su']:
            files = glob(self.cwd + '/' + 'traces/adj/*SU')
            unix.rename('_SU', '_SU.adj', files)

        # workaround for SPECFEM3D's requirement that all components exist,
        # even ones not in use
        unix.cd(self.cwd + '/' + 'traces/adj')
        for iproc in range(PAR.NPROC):
            for channel in ['x', 'y', 'z']:
                src = '%d_d%s_SU.adj' % (iproc, PAR.CHANNELS[0])
                dst = '%d_d%s_SU.adj' % (iproc, channel)
                if not exists(dst):
                    unix.cp(src, dst)

    def rename_data(self):
        """ Works around conflicting data filename conventions
        """
        if PAR.FORMAT in ['SU', 'su']:
            files = glob(self.cwd + '/' + 'traces/adj/*SU')
            unix.rename('_SU', '_SU.adj', files)

    def write_parameters(self):
        unix.cd(self.cwd)
        solvertools.write_parameters(vars(PAR))

    def write_receivers(self):
        unix.cd(self.cwd)
        key = 'use_existing_STATIONS'
        val = '.true.'
        setpar(key, val)
        _, h = preprocess.load('traces/obs')
        solvertools.write_receivers(h.nr, h.rx, h.rz)

    def write_sources(self):
        unix.cd(self.cwd)
        _, h = preprocess.load(dir='traces/obs')
        solvertools.write_sources(vars(PAR), h)

    # postprocessing wrapper overload
    def smooth(self, input_path='', output_path='',
               parameters=[], span_h=0., span_v=0.):
        """ Smooths kernels by convolving them with a Gaussian.  Wrapper over 
            xsmooth_sem utility. 
            smooth() in base.py has the incorrect command line call, specfem 
            requires that NPROC be specified
        """
        if not exists(input_path):
            raise Exception

        if not exists(output_path):
            unix.mkdir(output_path)

        # apply smoothing operator
        unix.cd(self.cwd)
        print 'smoothing parameters ', self.parameters
        for name in parameters or self.parameters:
            print 'smoothing', name
            solver_call = " ".join([
                    PATH.SPECFEM_BIN + '/' + 'xsmooth_sem',  # ./bin/xsmooth_sem
                    str(span_h),  # SIGMA_H
                    str(span_v),  # SIGMA_V
                    name + '_kernel',  # KERNEL_NAME
                    input_path + '/',  # INPUT_DIR
                    output_path + '/',  # OUTPUT_DIR
                    '.false'  # USE_GPU
                    ])
            call_solver(system.mpiexec(), solver_call)
        print ''

        # rename output files
        files = glob(output_path+'/*')
        unix.rename('_smooth', '', files)

    def combine_vol_data(self, output_path='', quantity=''):
        """
        This does not work
        Call Specfems executable combine_vol_data_vtk on kernels or model files
        """
        if not exists(output_path):
            unix.mkdir(output_path)
       
        # This should probably be moved to its own function 
        # def import_kernels()
        unix.cd(self.cwd)
        src = glob(join(PATH.GRAD, self.source_name, "*{}*".format(quantity)))
        dst = join(self.cwd, "kernels")
        unix.mkdir(dst)
        unix.ln(src=src, dst=dst)
        
        solver_call = " ".join([
                PATH.SPECFEM_BIN + '/' + 'xcombine_vol_data_vtk',
                0, # NPROC_START
                PAR.NPROC,  # NPROC_END
                quantity,  # QUANTITY
                dst,  # DIR_IN
                dst,  # DIR_OUT, we will rename the files first
                0  # GPU ACCEL
                ])
        call_solver(system_mpiexec(), solver_call)
        
        unix.rm(dst)
        print ''

    # miscellaneous
    @property
    def data_wildcard(self):
        channels = PAR.CHANNELS
        return '*_d[%s]_SU' % channels.lower()

    @property
    def data_filenames(self):
        unix.cd(self.cwd+'/'+'traces/obs')

        if PAR.FORMAT in ['SU', 'su']:
            if not PAR.CHANNELS:
                return sorted(glob('*_d?_SU'))
            filenames = []
            for channel in PAR.CHANNELS:
                filenames += sorted(glob('*_d'+channel+'_SU'))
            return filenames

        else:
            raise NotImplementedError

    @property
    def kernel_databases(self):
        return join(self.cwd, 'OUTPUT_FILES/DATABASES_MPI')

    @property
    def model_databases(self):
        return join(self.cwd, 'OUTPUT_FILES/DATABASES_MPI')

    @property
    def source_prefix(self):
        return 'CMTSOLUTION'

