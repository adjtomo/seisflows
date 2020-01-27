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

from glob import glob
from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import custom_import
from seisflows.tools.seismic import call_solver, getpar, setpar

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']


class Specfem3DPyatoa(custom_import('solver', 'base')):
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
        super(Specfem3DPyatoa, self).check()

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
            print(f"{self.__class__.__name__}.generate_data()")

        unix.cd(self.cwd)
        setpar("SIMULATION_TYPE", "1")
        setpar("SAVE_FORWARD", ".true.")
        setpar("ATTENUATION ", ".true.")
        call_solver(mpiexec=system.mpiexec(), executable="bin/xspecfem3D")

        # ASCII .sem? output format
        if PAR.FORMAT == "ascii":
            unix.mv(src=glob(os.path.join("OUTPUT_FILES", "*sem?")),
                    dst="traces/obs")
        else:
            raise NotImplementedError("Pyatoa only works with ascii")

        # Export traces to permanent storage on disk
        if PAR.SAVETRACES:
            self.export_traces(os.path.join(PATH.OUTPUT, "traces", "obs"))

    def generate_mesh(self, model_path, model_name, model_type='gll', 
                      symlink=True):
        """
        Performs meshing and database generation

        :type model_path: str
        :param model_path: path to the model to be used for mesh generation
        :type model_name: str
        :param model_name: name of the model to be used as identification
        :type model_type: str
        :param model_type: available model types to be passed to the Specfem3D
            Par_file. See Specfem3D Par_file for available options.
        :type symlink: bool
        :param symlink: symlink critical components rather than copying.
            This saves on storage requirements
        """
        # Determine which function to use when copying from Specfem directory
        if symlink:
            copy_func = unix.ln
        else:
            copy_func = unix.cp

        available_model_types = ["gll"]

        if PAR.VERBOSE:
            print(f"{self.__class__.__name__}.generate_mesh()")

        unix.cd(self.cwd)

        # Run mesh generation
        assert(exists(model_path))
        if model_type in available_model_types:
            par = getpar("MODEL").strip()
            if par == "gll":
                self.check_mesh_properties(model_path)

                src = glob(os.path.join(model_path, "*"))
                dst = self.model_databases
                copy_func(src, dst)

                call_solver(mpiexec=system.mpiexec(),
                            executable="bin/xgenerate_databases")
                
                # Remove VT? files, they aren't necessary and take up space
                for fid in glob(os.path.join(self.model_databases, "*.vt?")):
                    os.remove(fid)

            # Export the model for future use in the workflow
            if self.taskid == 0:
                self.export_model(os.path.join(PATH.OUTPUT, model_name))
        else:
            raise NotImplementedError(f"MODEL={model_type} not implemented")

    def eval_fwd(self, path=''):
        """
        Performs forward simulations for misfit function evaluation.
        Same as solver.base.eval_func without the residual writing.

        Function evaluation is taken care of by Pyatoa.
        """
        if PAR.VERBOSE:
            print(f"{self.__class__.__name__}.eval_fwd()")

        unix.cd(self.cwd)
        self.import_model(path)
        self.forward()
    
    def eval_func(self, pyaflowa):
        """
        Call Pyatoa workflow to evaluate the misfit functional.

        :type pyaflowa: Pyaflowa object
        :param pyaflowa: Pyaflowa object that controls the misfit function eval
        """
        if PAR.VERBOSE:
            print(f"{self.__class__.__name__}.eval_func() => calling Pyatoa...")
        pyaflowa.process(cwd=self.cwd, event_id=self.source_name)
    
    def forward(self, path='traces/syn'):
        """
        Overwrites seisflows.solver.specfem3d.forward()

        Calls SPECFEM3D forward solver and then moves files into path
        Adds attenuation, and output format for ASCII

        :type path: str
        :param path: path to export traces to after completion of simulation
        """
        # Set parameters and run forward simulation
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        setpar('ATTENUATION ', '.true.')
        call_solver(mpiexec=system.mpiexec(),
                    executable='bin/xgenerate_databases')
        call_solver(mpiexec=system.mpiexec(), executable='bin/xspecfem3D')

        # Seismic unix output format
        if PAR.FORMAT in ['SU', 'su']:
            unix.mv(src=glob(os.path.join("OUTPUT_FILES", "*_d?_SU")), dst=path)

        # ascii sem output format
        elif PAR.FORMAT == "ascii":
            unix.mv(src=glob(os.path.join("OUTPUT_FILES", "*sem?")), dst=path)

    def check_solver_parameter_files(self):
        """
        Checks solver parameters
        """
        nt = getpar(key="NSTEP", cast=int)
        dt = getpar(key="DT", cast=float)

        if nt != PAR.NT:
            if self.taskid == 0:
                warnings.warn("Specfem3D NSTEP != PAR.NT\n"
                              "overwriting Specfem3D with Seisflows parameter"
                              )
            setpar(key="NSTEP", val=PAR.NT)

        if dt != PAR.DT:
            if self.taskid == 0:
                warnings.warn("Specfem3D DT != PAR.DT\n"
                              "overwriting Specfem3D with Seisflows parameter"
                              )
            setpar(key="DT", val=PAR.DT)

        if self.mesh_properties.nproc != PAR.NPROC:
            if self.taskid == 0:
                warnings.warn("Specfem3D mesh nproc != PAR.NPROC")

        if 'MULTIPLES' in PAR:
            raise NotImplementedError

    def combine_vol_data_vtk(self):
        """
        Overwrites seisflows.solver.base.combine_vol_data_vtk

        Postprocessing wrapper: xcombine_vol_data_vtk

        Note:
            Default output of xcombine_vol_data_vtk will be named {quantity}.vtk
            This will create vtk files for all event kernels

        !!! should this be moved into base?
        """
        # Create a directory for combining kernels
        unix.cd(self.cwd)
        sum_dir = os.path.join(self.cwd, "SUM")
        if not os.path.exists(sum_dir):
            os.makedirs(sum_dir)

        # Reused solver call
        solver_call = " ".join([
            f"{PATH.SPECFEM_BIN}/xcombine_vol_data_vtk",
            f"0",  # NPROC_START
            f"{PAR.NPROC}",  # NPROC_END
            "{kernel}",  # QUANTITY
            "{DIR_IN}",  # DIR_IN
            "{DIR_OUT}",  # DIR_OUT
            f"0"  # GPU ACCEL
        ])

        # Parameters file determines which kernels are made into VTK files
        if PAR.VTK_EVENT_KERNELS:
            for kernel in PAR.VTK_EVENT_KERNELS:
                # Symlink kernels into the output path
                events = glob(os.path.join(PATH.GRAD, "kernels", "*"))
                for event in events:
                    # Deal with summed kernels separately
                    if event not in ["sum", "sum_nosmooth"]:
                        # symlink files into the SUM directory
                        unix.ln(src=glob(os.path.join(event, f"*{kernel}*")),
                                dst=sum_dir)
                        # Run the xcombine_vol_data_vtk
                        call_solver(mpiexec=system.mpiexec(),
                                    executable=solver_call.format(kernel=kernel,
                                                                  dir_in=event,
                                                                  dir_out=event)
                                    )
                        # Unfinished, how do I remove the symlinks and leave the
                        # VTK file, I need to rename it also and move it away
                        # Can I wrap this into a function?
                        unix.rm(glob(os.path.join(event, f"*{kernel}*")))

    @property
    def kernel_databases(self):
        """
        The location of databases for mernel outputs
        """
        return os.path.join(self.cwd, "OUTPUT_FILES", "DATABASES_MPI")

    @property
    def model_databases(self):
        """
        The location of databases for model outputs
        """
        return os.path.join(self.cwd, "OUTPUT_FILES", "DATABASES_MPI")

    @property  
    def source_prefix(self):
        """
        Specfem3D's preferred source prefix
        
        :rtype: str
        :return: source prefix
        """
        return "CMTSOLUTION"
       
          
