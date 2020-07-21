#!/usr/bin/env python
"""
This is the subclass seisflows.solver.Specfem3DSerial

This class is almost the same as Specfem3D, except the setup step is run as
a serial task. This is useful if HPC job queues are long, as it saves
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


class Specfem3DSerial(custom_import("solver", "specfem3d")):
    """
    Python interface to Specfem3D Cartesian. This subclass inherits functions
    from seisflows.solver.specfem3d

    !!! See base class for method descriptions !!!
    """
    def setup(self, model="init"):
        """
        Overload of solver.base.setup(), should be run as a single instance

        :type model: str
        :param model: "init", "true", generates the mesh to be used for workflow
            "true" used for synthetic-synthetic cases
            "init" for initial model, default
        :type model: str
        :param model: model to setup, either 'true' or 'init'
        """   
        if model == "true":
            self.generate_mesh(model_path=PATH.MODEL_TRUE,
                               model_name="model_true", model_type="gll")
        elif model == "init":
            self.generate_mesh(model_path=PATH.MODEL_INIT, 
                               model_name="model_init", model_type="gll")
            
        self.distribute_databases()

    def generate_data(self, **model_kwargs):
        """
        Overload seisflows.solver.base.generate_data. To be run in parallel
        
        Not used if PAR.CASE == "Data"

        Generates data in the synthetic-synthetic comparison case.
        Automatically calls generate mesh for the true model, rather than
        passing them in as kwargs.

        Also turns on attenuation for the forward model
        !!! attenuation could be moved into parameters.yaml? !!!
        """
        unix.cd(self.cwd)

        setpar("SIMULATION_TYPE", "1")
        setpar("SAVE_FORWARD", ".false.")
        setpar("ATTENUATION ", ".true.")

        call_solver(mpiexec=system.mpiexec(), executable="bin/xspecfem3D")

        # move ASCII .sem? files into appropriate directory
        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard)),
                dst=os.path.join("traces", "obs"))

        # Export traces to permanent storage on disk
        if PAR.SAVETRACES:
            self.export_traces(os.path.join(PATH.OUTPUT, "traces", "obs"))

    def generate_mesh(self, model_path, model_name, model_type='gll'):
        """
        Performs meshing and database generation as a serial task. Differs
        slightly from specfem3d class as it only creates database files for
        the main solver, rather than generating them individually for each
        running process.

        Note:
            Specfem model files must exist in each individual solver directory
            If symlinked, Fortran I/O errors will occur when multiple processes 
            try to read from the same .bin file. 
            This means that during the workflow there will be large disk 
            requirements due to redundant files in scratch.

        :type model_path: str
        :param model_path: path to the model to be used for mesh generation
        :type model_name: str
        :param model_name: name of the model to be used as identification
        :type model_type: str
        :param model_type: available model types to be passed to the Specfem3D
            Par_file. See Specfem3D Par_file for available options.
        """
        assert(exists(model_path)), f"model {model_path} does not exist"

        available_model_types = ["gll"]
        assert(model_type in available_model_types), \
            f"{model_type} not in available types {available_model_types}"

        # Ensure this is only run by taskid 0, mainsolver
        cwd = os.path.join(PATH.SOLVER, self.mainsolver)
        unix.cd(cwd)

        par = getpar("MODEL").strip()
        if par == "gll":
            self.check_mesh_properties(model_path)
            
            # Copy model files and then run xgenerate databases
            src = glob(os.path.join(model_path, "*"))
            dst = self.model_databases
            unix.cp(src, dst)
            call_solver(mpiexec=system.mpiexec(),
                        executable="bin/xgenerate_databases")

        self.export_model(os.path.join(PATH.OUTPUT, model_name))

    def distribute_databases(self):
        """
        A serial task to distrubute the database files outputted by 
        xgenerate_databases from main solver to all solver directories
        """
        # Copy the database files but ignore any vt? files
        src_db = glob(os.path.join(PATH.SOLVER, self.mainsolver, 
                                   "OUTPUT_FILES", "DATABASES_MPI", "*"))
        for extension in [".vtu", ".vtk"]:
            src_db = [_ for _ in src_db if extension not in _] 
    
        # Copy the .h files from the mesher, Specfem needs these as well
        src_h = glob(os.path.join(PATH.SOLVER, self.mainsolver, 
                                  "OUTPUT_FILES", "*.h"))

        for source_name in self.source_names:
            # Ensure main solver is skipped
            if source_name == self.mainsolver:
                continue
            # Copy database files
            dst_db = os.path.join(PATH.SOLVER, source_name, 
                                  "OUTPUT_FILES", "DATABASES_MPI", "")
            unix.cp(src_db, dst_db)

            # Copy mesher h files
            dst_h = os.path.join(PATH.SOLVER, source_name, "OUTPUT_FILES", "")
            unix.cp(src_h, dst_h)

    def initialize_solver_directories(self):
        """
        Overwrite solver.base.initialize_solver_directories()

        Creates solver directories in serial.

        Should only be run by master job, copies the necessary data and creates
        the necessary file structure for all events.
        """
        for source_name in self.source_names:
            cwd = os.path.join(PATH.SOLVER, source_name)
            # Remove any existing scratch directory 
            unix.rm(cwd)

            # Create internal directory structure, change into directory to make
            # all actions RELATIVE path actions
            unix.mkdir(cwd)
            unix.cd(cwd)
            for cwd_dir in ["bin", "DATA", "OUTPUT_FILES/DATABASES_MPI", 
                            "traces/obs", "traces/syn", "traces/adj"]:
                unix.mkdir(cwd_dir)

            # Copy exectuables
            src = glob(os.path.join(PATH.SPECFEM_BIN, "*"))
            dst = os.path.join("bin", "")
            unix.cp(src, dst)

            # Copy all input files except source files
            src = glob(os.path.join(PATH.SPECFEM_DATA, "*"))
            src = [_ for _ in src if not self.source_prefix in _]
            dst = os.path.join("DATA", "")
            unix.cp(src, dst)

            # symlink event source specifically
            src = os.path.join(PATH.SPECFEM_DATA, 
                               f"{self.source_prefix}_{source_name}")
            dst = os.path.join("DATA", self.source_prefix)
            unix.ln(src, dst)

            if source_name == self.mainsolver: 
                # Symlink taskid_0 as mainsolver in solver directory
                unix.ln(source_name, os.path.join(PATH.SOLVER, "mainsolver"))
                # Only check the solver parameters once
                self.check_solver_parameter_files()

    @property
    def mainsolver(self):
        """
        Ensure that the main solver has a consistent reference inside Solver
        """
        return self.source_names[0]

