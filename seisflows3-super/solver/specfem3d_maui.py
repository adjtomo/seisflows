#!/usr/bin/env python3
"""
This is the subclass seisflows.solver.Specfem3DMaui

This class is almost the same as Specfem3D, except the setup step is run as
a serial task. This is useful as HPC job queues are long on Maui, so it saves
on job queue time by replacing it with a long-winded serial task.

Additionally, misfit quantification is split off from the forward simulation,
because Anaconda is not available on the main cluster, so jobs need to be
submitted to an auxiliary cluster. This is paired with a new evaluate_function()
function defined by the InversionMaui workflow class.
"""
import os
import sys
import warnings

from glob import glob
from seisflows3.tools import unix
from seisflows3.tools.wrappers import exists
from seisflows3.config import custom_import, SeisFlowsPathsParameters
from seisflows3.tools.specfem import call_solver, getpar, setpar


# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']


class Specfem3DMaui(custom_import("solver", "specfem3d")):
    """
    Python interface to Specfem3D Cartesian. This subclass inherits functions
    from seisflows3.solver.specfem3d

    !!! See base class for method descriptions !!!
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

    def setup(self, model):
        """
        Overload of solver.base.setup(), should be run as a single instance

        :type model: str
        :param model: "init", "true", generates the mesh to be used for workflow
            "true" used for synthetic-synthetic cases
            "init" for initial model, default
        :type model: str
        :param model: model to setup, either 'true' or 'init'
        """
        # Choice of model will determine which mesh to generate
        self.generate_mesh(model_path=getattr(PATH, f"MODEL_{model.upper()}"),
                           model_name=f"model_{model.lower()}",
                           model_type="gll")
            
        self.distribute_databases()

    def generate_data(self):
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

        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")
        if PAR.ATTENUATION:
            setpar(key="ATTENUATION ", val=".true.", file="DATA/Par_file")
        else:
            setpar(key="ATTENUATION ", val=".false.", file="DATA/Par_file")

        call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xspecfem3D")

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
        the main solver, which are then copied in serial by the function
        distribute_databases()

        :type model_path: str
        :param model_path: path to the model to be used for mesh generation
        :type model_name: str
        :param model_name: name of the model to be used as identification
        :type model_type: str
        :param model_type: available model types to be passed to the Specfem3D
            Par_file. See Specfem3D Par_file for available options.
        """
        available_model_types = ["gll"]

        assert(exists(model_path)), f"model {model_path} does not exist"

        model_type = model_type or getpar(key="MODEL", file="DATA/Par_file")
        assert(model_type in available_model_types), \
            f"{model_type} not in available types {available_model_types}"

        # Ensure that we're running on the main solver only
        assert(self.taskid == 0)

        unix.cd(self.cwd)

        # Check that the model parameter falls into the acceptable types
        par = getpar("MODEL").strip()
        assert(par in available_model_types), \
            f"Par_file {par} not in available types {available_model_types}"

        if par == "gll":
            self.check_mesh_properties(model_path)
            
            # Copy model files and then run xgenerate databases
            src = glob(os.path.join(model_path, "*"))
            dst = self.model_databases
            unix.cp(src, dst)

            call_solver(mpiexec=PAR.MPIEXEC,
                        executable="bin/xgenerate_databases")

        self.export_model(os.path.join(PATH.OUTPUT, model_name))

    def eval_misfit(self, path='', export_traces=False):
        """
        Performs function evaluation only, that is, the misfit quantification.
        Forward simulations are performed in a separate function

        :type path: str
        :param path: path in the scratch directory to use for I/O
        :type export_traces: bool
        :param export_traces: option to save the observation traces to disk
        :return:
        """
        preprocess.prepare_eval_grad(cwd=self.cwd, taskid=self.taskid,
                                     source_name=self.source_name,
                                     filenames=self.data_filenames
                                     )
        if export_traces:
            self.export_residuals(path)

    def eval_fwd(self, path=''):
        """
        High level solver interface

        Performans forward simulations only, function evaluation is split off
        into its own function

        :type path: str
        :param path: path in the scratch directory to use for I/O
        """
        unix.cd(self.cwd)
        self.import_model(path)
        self.forward()

    def distribute_databases(self):
        """
        A serial task to distrubute the database files outputted by 
        xgenerate_databases from main solver to all other solver directories
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
            # Copy database files to each of the other source directories
            dst_db = os.path.join(PATH.SOLVER, source_name, 
                                  "OUTPUT_FILES", "DATABASES_MPI", "")
            unix.cp(src_db, dst_db)

            # Copy mesher h files into the overlying directory
            dst_h = os.path.join(PATH.SOLVER, source_name, "OUTPUT_FILES", "")
            unix.cp(src_h, dst_h)

    def initialize_solver_directories(self):
        """
        Creates solver directories in serial using a single node.
        Should only be run by master job.

        Differs from Base initialize_solver_directories() as this serial task
        will create directory structures for each source, rather than having
        each source create its own. However the internal dir structure is the
        same.
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
            src = [_ for _ in src if self.source_prefix not in _]
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

    def check_solver_parameter_files(self):
        """
        Checks solver parameters. Only slightly different to Specfem3D as it
        is run by the main task, not be an array process, so no need to check
        task_id
        """
        nt = getpar(key="NSTEP", cast=int)
        dt = getpar(key="DT", cast=float)

        if nt != PAR.NT:
            warnings.warn("Specfem3D NSTEP != PAR.NT\n"
                          "overwriting Specfem3D with Seisflows parameter"
                          )
            setpar(key="NSTEP", val=PAR.NT)

        if dt != PAR.DT:
            warnings.warn("Specfem3D DT != PAR.DT\n"
                          "overwriting Specfem3D with Seisflows parameter"
                          )
            setpar(key="DT", val=PAR.DT)

        if self.mesh_properties.nproc != PAR.NPROC:
            warnings.warn("Specfem3D mesh nproc != PAR.NPROC")

        if "MULTIPLES" in PAR:
            raise NotImplementedError

    @property
    def mainsolver(self):
        """
        Ensure that the main solver has a consistent reference inside Solver
        """
        return self.source_names[0]

