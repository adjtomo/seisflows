#!/usr/bin/env python3
"""
This Solver module is in charge of interacting with external numerical solvers
such as SPECFEM (2D/3D/3D_GLOBE). The Base class provides general functions
that work with SPECFEM, while subclasses provide details to differentiate the
various types of SPECFEM.
"""
import os
import sys
import logging
import subprocess
import numpy as np
from glob import glob
from functools import partial

from seisflows.plugins import solver_io
from seisflows.tools import msg, unix
from seisflows.tools.specfem import Container
from seisflows.tools.wrappers import Struct, diff, exists
from seisflows.config import SeisFlowsPathsParameters


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']
system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']


class Base:
    """
    This base class provides an interface through which solver simulations can
    be set up and run and a parent class for the following subclasses:

    SPECFEM2D
    SPECFEM3D
    SPECFEM3D_GLOBE

    .. note:::
        This class supports only acoustic and isotropic elastic inversions.

    Function descriptors:

    eval_func, eval_grad, apply_hess

        These methods deal with evaluation of the misfit function or its
        derivatives.  Together, they provide the primary interface through which
        SeisFlows interacts with SPECFEM2D/3D

    forward, adjoint

        These methods allow direct access to low-level SPECFEM2D/3D components,
        providing an alternative interface through which to interact with the
        solver

    steup, generate_data, generate_model

        One-time operations performed at the beginning of inversion or migration

    initialize_solver_directories, initialize_adjoint_traces

        SPECFEM2D/3D requires a particular directory structure in which to run
        and particular file formats for models, data, and parameter files. These
        methods help put in place all these prerequisites

    load, save

        For reading and writing SPECFEM2D/3D models and kernels. On the disk,
        models and kernels are stored as binary files, and in memory, as
        dictionaries with different keys corresponding to different material
        parameters

    split, merge

        Within the solver routines, it is natural to store models as
        dictionaries. Within the optimization routines, it is natural to store
        models as vectors. Two methods, 'split' and 'merge', are used to convert
        back and forth between these two representations

    combine, smooth

        Utilities for combining and smoothing kernels

    generate_data, generate_mesh, eval_fwd, forward, adjoint,
    data_filenames, model_databases, kernel_databases, source_prefix

        !!! Required functions which must be implemented by subclass !!!

    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the User_!
        Attributes are just initialized as NoneTypes for clarity and docstrings

        :type parameters: list of str
        :param parameters: a list detailing the parameters to be used to
            define the model, available: ['vp', 'vs', 'rho']
        :type _mesh_properties: seisflows.tools.wrappers.Struct
        :param _mesh_properties: hidden attribute, a dictionary of mesh
            properties, including the ngll points, nprocs, and mesh coordinates
        :type _source_names: hidden attribute,
        :param _source_names: the names of all the sources that are being used
            by the solver
        :type logger: Logger
        :param logger: Class-specific logging module, log statements pushed
            from this logger will be tagged by its specific module/classname
        """
        self.parameters = []
        self._mesh_properties = None
        self._source_names = None

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        sf.par("MATERIALS", required=True, par_type=str,
               docstr="Material parameters used to define model. Available: "
                      "['ELASTIC': Vp, Vs, 'ACOUSTIC': Vp, 'ISOTROPIC', "
                      "'ANISOTROPIC']")

        sf.par("DENSITY", required=True, par_type=str,
               docstr="How to treat density during inversion. Available: "
                      "['CONSTANT': Do not update density, "
                      "'VARIABLE': Update density]")

        sf.par("ATTENUATION", required=True, par_type=bool,
               docstr="If True, turn on attenuation during forward "
                      "simulations, otherwise set attenuation off. Attenuation "
                      "is always off for adjoint simulations.")

        sf.par("COMPONENTS", required=False, default="ZNE", par_type=str,
               docstr="Components used to generate data, formatted as a single "
                      "string, e.g. ZNE or NZ or E")

        sf.par("SOLVERIO", required=False, default="fortran_binary",
               par_type=int,
               docstr="The format external solver files. Available: "
                      "['fortran_binary', 'adios']")

        sf.path("SOLVER", required=False,
                default=os.path.join(PATH.SCRATCH, "solver"),
                docstr="scratch path to hold solver working directories")

        sf.path("SPECFEM_BIN", required=True,
                docstr="path to the SPECFEM binary executables")

        sf.path("SPECFEM_DATA", required=True,
                docstr="path to the SPECFEM DATA/ directory containing the "
                       "'Par_file', 'STATIONS' file and 'CMTSOLUTION' files")

        sf.path("DATA", required=False, 
                docstr="path to a directory containing any external data "
                       "required by the workflow. Catch all directory that "
                       "can be accessed by all modules")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        if validate:
            self.required.validate()

        # Check that other modules have set parameters that will be used here
        for required_parameter in ["NPROC"]:
            assert (required_parameter in PAR), \
                f"Solver requires {required_parameter}"

        available_materials = ["ELASTIC", "ACOUSTIC",  # specfem2d, specfem3d
                               "ISOTROPIC", "ANISOTROPIC"]  # specfem3d_globe
        assert(PAR.MATERIALS.upper() in available_materials), \
            f"MATERIALS must be in {available_materials}"

        acceptable_densities = ["CONSTANT", "VARIABLE"]
        assert(PAR.DENSITY.upper() in acceptable_densities), \
            f"DENSITY must be in {acceptable_densities}"

        # Internal parameter list based on user-input material choices
        # Important to reset parameters to a blank list and let the check
        # statements fill it. If not, each time workflow is resumed, parameters
        # list will append redundant parameters and things stop working
        self.parameters = []
        if PAR.MATERIALS.upper() == "ELASTIC":
            assert(PAR.SOLVER.lower() in ["specfem2d", "specfem3d"])
            self.parameters += ["vp", "vs"]
        elif PAR.MATERIALS.upper() == "ACOUSTIC":
            assert(PAR.SOLVER.lower() in ["specfem2d", "specfem3d"])
            self.parameters += ["vp"]
        elif PAR.MATERIALS.upper() == "ISOTROPIC":
            assert(PAR.SOLVER.lower() in ["specfem3d_globe"])
            self.parameters += ["vp", "vs"]
        elif PAR.MATERIALS.upper() == "ANISOTROPIC":
            assert(PAR.SOLVER.lower() in ["specfem3d_globe"])
            self.parameters += ["vpv", "vph", "vsv", "vsh", "eta"]

        if PAR.DENSITY.upper() == "VARIABLE":
            self.parameters.append("rho")

        assert hasattr(solver_io, PAR.SOLVERIO)
        assert hasattr(self.io, "read_slice"), \
            "IO method has no attribute 'read_slice'"
        assert hasattr(self.io, "write_slice"), \
            "IO method has no attribute 'write_slice'"

    def setup(self):
        """ 
        Prepares solver for inversion or migration.
        Sets up directory structure expected by SPECFEM and copies or generates
        seismic data to be inverted or migrated

        .. note:;
            As input for an inversion or migration, users can choose between
            providing data, or providing a target model from which data are
            generated on the fly.
            In the former case, a value for PATH.DATA must be supplied;
            in the latter case, a value for PATH.MODEL_TRUE must be provided.
        """
        unix.rm(self.cwd)
        self.initialize_solver_directories()
        self.check_solver_parameter_files()
        self.generate_data()
        self.generate_mesh(model_name="init", model_type="gll")
        self.initialize_adjoint_traces()

    def clean(self):
        """
        Clean up solver-dependent run directory by removing the OUTPUT_FILES/
        directory
        """
        unix.cd(self.cwd)
        unix.rm("OUTPUT_FILES")
        unix.mkdir("OUTPUT_FILES")

    def generate_mesh(self, model_path, model_name, model_type):
        """
        Performs meshing and database generation.

        This function is Solver specific and is responsible for generating
        the mesh using the external numerical solver.
        """
        raise NotImplementedError("must be implemented by solver subclass")

    def generate_data(self):
        """
        Performs meshing and data generation for "true" data.
        """
        raise NotImplementedError("must be implemented by solver subclass")

    def eval_func(self, path, write_residuals=True):
        """
        Performs forward simulations and evaluates the misfit function

        :type path: str
        :param path: directory from which model is imported and where residuals
            will be exported
        :type write_residuals: bool
        :param write_residuals: calculate and export residuals
        """
        if self.taskid == 0:
            self.logger.info("running forward simulations")

        unix.cd(self.cwd)
        self.import_model(path)
        self.forward()

        if write_residuals:
            if self.taskid == 0:
                self.logger.debug("calling preprocess.prepare_eval_grad()")
            preprocess.prepare_eval_grad(cwd=self.cwd, taskid=self.taskid,
                                         source_name=self.source_name,
                                         filenames=self.data_filenames
                                         )
            self.export_residuals(path)

    def eval_grad(self, path, export_traces=False):
        """
        High level solver interface that evaluates gradient by carrying out
        adjoint simulations. A function evaluation must already have been
        carried out.

        :type path: str
        :param path: directory from which model is imported
        :type export_traces: bool
        :param export_traces: if True, save traces to OUTPUT.
            if False, discard traces
        """
        unix.cd(self.cwd)
        if self.taskid == 0:
            self.logger.debug("running adjoint simulations")

        # Check to make sure that preprocessing module created adjoint traces
        adjoint_traces_wildcard = os.path.join("traces", "adj", "*")
        if not glob(adjoint_traces_wildcard):
            print(msg.cli(f"Event {self.source_name} has no adjoint traces, "
                          f"which will lead to an external solver error. "
                          f"Please check that solver.eval_func() executed "
                          f"properly", border="=", header="solver error")
                  )
            sys.exit(-1)

        self.adjoint()
        self.export_kernels(path)

        if export_traces:
            self.export_traces(path=os.path.join(path, "traces", "syn"),
                               prefix="traces/syn")
            self.export_traces(path=os.path.join(path, "traces", "adj"),
                               prefix="traces/adj")

    def apply_hess(self, path):
        """
        High level solver interface that computes action of Hessian on a given
        model vector. A gradient evaluation must have already been carried out.

        TODO preprocess has no function prepare_apply_hess()

        :type path: str
        :param path: directory to which output files are exported
        """
        raise NotImplementedError("must be implemented by solver subclass")

        unix.cd(self.cwd)
        self.import_model(path)
        unix.mkdir("traces/lcg")

        self.forward("traces/lcg")
        preprocess.prepare_apply_hess(self.cwd)
        self.adjoint()
        self.export_kernels(path)

    def forward(self, path):
        """
        Low level solver interface

        Calls forward solver

        !!! Must be implemented by subclass !!!
        """
        raise NotImplementedError("must be implemented by solver subclass")

    def adjoint(self):
        """
        Low level solver interface

        Calls adjoint solver

        !!! Must be implemented by subclass !!!
        """
        raise NotImplementedError("must be implemented by solver subclass")

    def call_solver(self, executable, output="solver.log"):
        """
        Calls MPI solver executable to run solver binaries, used by individual
        processes to run the solver on system. If the external solver returns a
        non-zero exit code (failure), this function will return a negative boolean.

        :type mpiexec: str
        :param mpiexec: call to mpi. If None (e.g., serial run, defaults to ./)
        :type executable: str
        :param executable: executable function to call
        :type output: str
        :param output: where to redirect stdout
        """
        if not os.path.exists(executable):
            print(msg.cli(f"solver executable {executable} does not exist",
                          header="external solver error", border="="))
            sys.exit(-1)

        # mpiexec is None when running in serial mode, so e.g., ./xmeshfem2D
        if PAR.SYSTEM in ["workstation"]:
            exc_cmd = f"./{executable}"
        # Otherwise mpiexec is system dependent (e.g., srun, mpirun)
        else:
            exc_cmd = f"{PAR.MPIEXEC} {executable}"

        try:
            # Write solver stdout (log files) to text file
            f = open(output, "w")
            subprocess.run(exc_cmd, shell=True, check=True, stdout=f)
        except (subprocess.CalledProcessError, OSError) as e:
            print(msg.cli("The external numerical solver has returned a nonzero "
                          "exit code (failure). Consider stopping any currently "
                          "running jobs to avoid wasted computational resources. "
                          f"Check 'scratch/solver/mainsolver/{output}' for the "
                          f"solvers stdout log message. "
                          f"The failing command and error message are: ",
                          items=[f"exc: {exc_cmd}", f"err: {e}"],
                          header="external solver error",
                          border="=")
                  )
            sys.exit(-1)
        finally:
            f.close()

    @property
    def io(self):
        """
        Solver IO module set by User.

        Located in seisflows.plugins.solver_io
        """
        return getattr(solver_io, PAR.SOLVERIO)

    def load(self, path, prefix="", suffix="", parameters=None,):
        """ 
        Solver I/O: Loads SPECFEM2D/3D models or kernels

        :type path: str
        :param path: directory from which model is read
        :type prefix: str
        :param prefix: optional filename prefix
        :type suffix: str
        :param suffix: optional filename suffix, eg '_kernel'
        :type parameters: list
        :param parameters: material parameters to be read
            (if empty, defaults to self.parameters)
        :rtype: dict
        :return: model or kernels indexed by material parameter and
            processor rank, ie dict[parameter][iproc]
        """
        if parameters is None:
            parameters = self.parameters

        load_dict = Container()
        for iproc in range(self.mesh_properties.nproc):
            for key in parameters:
                load_dict[key] += self.io.read_slice(
                    path=path, parameters=f"{prefix}{key}{suffix}", iproc=iproc
                )

        return load_dict

    def save(self, save_dict, path, parameters=None, prefix="", suffix=""):
        """ 
        Solver I/O: Saves SPECFEM2D/3D models or kernels

        :type save_dict: dict or Container
        :param save_dict: model stored as a dictionary or Container
        :type path: str
        :param path: directory from which model is read
        :type parameters: list
        :param parameters: list of material parameters to be read
        :type prefix: str
        :param prefix: optional filename prefix
        :type suffix: str
        :param suffix: optional filename suffix, eg '_kernel'
        """
        unix.mkdir(path)

        if parameters is None:
            parameters = self.parameters

        # Fill in any missing parameters
        missing_keys = diff(parameters, save_dict.keys())
        for iproc in range(self.mesh_properties.nproc):
            for key in missing_keys:
                save_dict[key] += self.io.read_slice(
                    path=PATH.MODEL_INIT, parameters=f"{prefix}{key}{suffix}",
                    iproc=iproc
                )

        # Write slices to disk
        for iproc in range(self.mesh_properties.nproc):
            for key in parameters:
                self.io.write_slice(data=save_dict[key][iproc], path=path,
                                    parameters=f"{prefix}{key}{suffix}",
                                    iproc=iproc)

    def merge(self, model, parameters=None):
        """
        Convert dictionary representation `model` to vector representation `m`

        :type model: dict
        :param model: model to be converted
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        :rtype: np.ndarray
        :return: model as a vector
        """
        if parameters is None:
            parameters = self.parameters

        m = np.array([])
        for key in parameters:
            for iproc in range(self.mesh_properties.nproc):
                m = np.append(m, model[key][iproc])

        return m

    def split(self, m, parameters=None):
        """
        Converts vector representation `m` to dictionary representation `model`

        :type m: np.ndarray
        :param m: model to be converted
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        :rtype: dict
        :return: model as a dictionary
        """
        if parameters is None:
            parameters = self.parameters

        nproc = self.mesh_properties.nproc
        ngll = self.mesh_properties.ngll
        model = Container()

        for idim, key in enumerate(parameters):
            model[key] = []
            for iproc in range(nproc):
                imin = sum(ngll) * idim + sum(ngll[:iproc])
                imax = sum(ngll) * idim + sum(ngll[:iproc + 1])
                model[key] += [m[imin:imax]]

        return model

    def combine(self, input_path, output_path, parameters=None):
        """
        Postprocessing wrapper: xcombine_sem 
        Sums kernels from individual source contributions to create gradient.

        .. note::
            The binary xcombine_sem simply sums matching databases (.bin) 

        .. note::
            It is ASSUMED that this function is being called by
            system.run(single=True) so that we can use the main solver
            directory to perform the kernel summation task

        :type input_path: str
        :param input_path: path to data
        :type output_path: str
        :param output_path: path to export the outputs of xcombine_sem
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        """
        if parameters is None:
            parameters = self.parameters

        if not exists(output_path):
            unix.mkdir(output_path)

        unix.cd(self.cwd)

        # Write the source names into the kernel paths file for SEM/ directory
        with open("kernel_paths", "w") as f:
            f.writelines(
                [os.path.join(input_path, f"{name}\n")
                 for name in self.source_names]
            )

        # Call on xcombine_sem to combine kernels into a single file
        for name in self.parameters:
            # e.g.: mpiexec ./bin/xcombine_sem alpha_kernel kernel_paths output
            self.call_solver(mpiexec=PAR.MPIEXEC,
                             executable=" ".join([
                                 f"bin/xcombine_sem", f"{name}_kernel",
                                 "kernel_paths", output_path]
                             )
                        )

    def smooth(self, input_path, output_path, parameters=None, span_h=0.,
               span_v=0., output="solver.log"):
        """
        Postprocessing wrapper: xsmooth_sem
        Smooths kernels by convolving them with a Gaussian.

        .. note::
            paths require a trailing `/` character when calling xsmooth_sem

        .. note::
            It is ASSUMED that this function is being called by
            system.run(single=True) so that we can use the main solver
            directory to perform the kernel smooth task

        :type input_path: str
        :param input_path: path to data
        :type output_path: str
        :param output_path: path to export the outputs of xcombine_sem
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        :type span_h: float
        :param span_h: horizontal smoothing length in meters
        :type span_v: float
        :param span_v: vertical smoothing length in meters
        :type output: str
        :param output: file to output stdout to
        """
        if parameters is None:
            parameters = self.parameters

        if not exists(output_path):
            unix.mkdir(output_path)

        # Apply smoothing operator inside scratch/solver/*
        unix.cd(self.cwd)

        # mpiexec ./bin/xsmooth_sem SMOOTH_H SMOOTH_V name input output use_gpu
        for name in parameters:
            self.call_solver(mpiexec=PAR.MPIEXEC,
                             executable=" ".join(["bin/xsmooth_sem",
                                                  str(span_h), str(span_v),
                                                  f"{name}_kernel",
                                                  os.path.join(input_path, ""),
                                                  os.path.join(output_path, ""),
                                                  ".false"]),
                             output=output
                        )

        # Rename output files
        files = glob(os.path.join(output_path, "*"))
        unix.rename(old="_smooth", new="", names=files)

    def import_model(self, path):
        """
        File transfer utility. Import the model into the workflow.

        :type path: str
        :param path: path to model
        """
        model = self.load(path=os.path.join(path, "model"))
        self.save(model, self.model_databases)

    def import_traces(self, path):
        """
        File transfer utility. Import traces into the workflow.

        :type path: str
        :param path: path to traces
        """
        src = glob(os.path.join(path, 'traces', self.source_name, '*'))
        dst = os.path.join(self.cwd, 'traces', 'obs')
        unix.cp(src, dst)

    def export_model(self, path, parameters=None):
        """
        File transfer utility. Export the model to disk.

        Performed by master solver.

        :type path: str
        :param path: path to save model
        :type parameters: list
        :param parameters: list of parameters that define the model
        """
        if parameters is None:
            parameters = self.parameters

        if self.taskid == 0:
            unix.mkdir(path)
            for key in parameters:
                files = glob(os.path.join(self.model_databases, f"*{key}.bin"))
                unix.cp(files, path)

    def export_kernels(self, path):
        """
        File transfer utility. Export kernels to disk

        :type path: str
        :param path: path to save kernels
        """
        if self.taskid == 0:
            self.logger.debug(f"exporting kernels to:\n{path}")

        unix.cd(self.kernel_databases)

        # Work around conflicting name conventions
        self.rename_kernels()

        src = glob("*_kernel.bin")
        dst = os.path.join(path, "kernels", self.source_name)
        unix.mkdir(dst)
        unix.mv(src, dst)

    def export_residuals(self, path):
        """
        File transfer utility. Export residuals to disk.

        :type path: str
        :param path: path to save residuals
        """
        if self.taskid == 0:
            self.logger.debug(f"exporting residuals to:\n{path}")

        unix.mkdir(os.path.join(path, "residuals"))
        src = os.path.join(self.cwd, "residuals")

        # If this residuals directory has not been created, something
        # has gone wrong with the preprocessing and workflow cannot proceed
        if not os.path.exists(src):
            print(msg.cli("The Solver function 'export_residuals' expected "
                          "'residuals' directories to be created but could not "
                          "find them and cannot continue the workflow. Please "
                          "check the preprocess.prepare_eval_grad() function",
                          header="preprocess error", border="="))
            sys.exit(-1)

        dst = os.path.join(path, "residuals", self.source_name)
        unix.mv(src, dst)

    def export_traces(self, path, prefix="traces/obs"):
        """
        File transfer utility. Export traces to disk.

        :type path: str
        :param path: path to save traces
        :type prefix: str
        :param prefix: location of traces w.r.t self.cwd
        """
        if self.taskid == 0:
            self.logger.debug("exporting traces to {path} {prefix}")

        unix.mkdir(os.path.join(path))

        src = os.path.join(self.cwd, prefix)
        dst = os.path.join(path, self.source_name)
        unix.cp(src, dst)

    def rename_kernels(self):
        """
        Works around conflicting kernel filename conventions by renaming
        `alpha` to `vp` and `beta` to `vs`
        """
        # Rename 'alpha' to 'vp'
        for tag in ["alpha", "alpha[hv]", "reg1_alpha", "reg1_alpha[hv]"]:
            names = glob(f"*proc??????_{tag}_kernel.bin")
            unix.rename(old="alpha", new="vp", names=names)

        # Rename 'beta' to 'vs'
        for tag in ["beta", "beta[hv]", "reg1_beta", "reg1_beta[hv]"]:
            names = glob(f"*proc??????_{tag}_kernel.bin")
            unix.rename(old="beta", new="vs", names=names)

    def rename_data(self, path):
        """
        Optional method to rename data to work around conflicting naming schemes
        for data outputted by the solver

        !!! Can be implemented by subclass !!!
        """
        raise NotImplementedError("must be implemented by solver subclass")

    def initialize_solver_directories(self):
        """
        Creates directory structure expected by SPECFEM3D (bin/, DATA/) copies
        executables, and prepares input files. Executables must be supplied
        by user as there is no mechanism for automatically compiling from source

        Directories will act as completely independent Specfem run directories.
        This allows for embarrassing parallelization while avoiding the need
        for intra-directory communications, at the cost of redundancy and
        extra files.
        """
        if self.taskid == 0:
            self.logger.info(f"initializing {PAR.NTASK} solver directories")

        unix.mkdir(self.cwd)
        unix.cd(self.cwd)

        # Create directory structure
        for cwd_dir in ["bin", "DATA", "OUTPUT_FILES/DATABASES_MPI",
                        "traces/obs", "traces/syn", "traces/adj",
                        self.model_databases, self.kernel_databases]:
            unix.mkdir(cwd_dir)

        # Copy exectuables into the bin/ directory
        src = glob(os.path.join(PATH.SPECFEM_BIN, "*"))
        dst = os.path.join("bin", "")
        unix.cp(src, dst)

        # Copy all input files except source files
        src = glob(os.path.join(PATH.SPECFEM_DATA, "*"))
        src = [_ for _ in src if self.source_prefix not in _]
        dst = os.path.join("DATA", "")
        unix.cp(src, dst)

        # Symlink event source specifically, strip the source name as SPECFEM
        # just expects `source_name`
        src = os.path.join(PATH.SPECFEM_DATA, 
                           f"{self.source_prefix}_{self.source_name}")
        dst = os.path.join("DATA", self.source_prefix)
        unix.ln(src, dst)

        if self.taskid == 0: 
            mainsolver = os.path.join(PATH.SOLVER, "mainsolver")
            # Symlink taskid_0 as mainsolver in solver directory for convenience
            if not os.path.exists(mainsolver):
                unix.ln(self.cwd, mainsolver)
                self.logger.debug(f"source {self.source_name} symlinked as "
                                  f"mainsolver")
        else:
            # Copy the initial model from mainsolver into current directory
            # Avoids the need to run multiple instances of xgenerate_databases
            src = os.path.join(PATH.SOLVER, "mainsolver", "OUTPUT_FILES", 
                               "DATABASES_MPI")
            dst = os.path.join(self.cwd, "OUTPUT_FILES", "DATABASES_MPI")
            unix.cp(src, dst)

    def initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM.
        This is only done for the 'base' the Preprocess class.

        .. note::
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        if PAR.PREPROCESS.lower() == "default":
            if self.taskid == 0:
                self.logger.debug(f"intializing {len(self.data_filenames)} "
                                  f"empty adjoint traces per event")

            for filename in self.data_filenames:
                st = preprocess.reader(
                            path=os.path.join(self.cwd, "traces", "obs"),
                            filename=filename
                            )
                # Zero out data just so we have empty adjoint traces as SPECFEM
                # will expect all adjoint sources to have all components
                st *= 0

                # Write traces back to the adjoint trace directory
                preprocess.writer(st=st, filename=filename,
                                  path=os.path.join(self.cwd, "traces", "adj")
                                  )

    def check_mesh_properties(self, path=None):
        """
        Determine if Mesh properties are okay for workflow

        TODO fix or rewrite this function

        :type path: str
        :param path: path to the mesh file
        """
        # Check the given model path or the initial model
        if path is None:
            path = PATH.MODEL_INIT

        if not exists(path):
            print(msg.cli(f"The following mesh path does not exist but should",
                          items=[path], header="solver error", border="="))
            sys.exit(-1)

        # Count slices and grid points
        key = self.parameters[0]
        iproc = 0
        ngll = []
        while True:
            dummy = self.io.read_slice(path=path, parameters=key, 
                                       iproc=iproc)[0]
            ngll += [len(dummy)]
            iproc += 1
            if not exists(os.path.join(path,
                                       f"proc{int(iproc):06d}_{key}.bin")):
                break
        nproc = iproc

        # Create coordinate pointers
        # !!! This partial is incorrectly defined and does not execute when 
        # !!! called. What is the point of that?
        coords = Struct()
        for key in ['x', 'y', 'z']:
            coords[key] = partial(self.io.read_slice, self, path, key)

        # Define internal mesh properties
        self._mesh_properties = Struct([["nproc", nproc],
                                        ["ngll", ngll],
                                        ["path", path],
                                        ["coords", coords]]
                                       )

    def check_source_names(self):
        """
        Determines names of sources by applying wildcard rule to
        user-supplied input files

        .. note::
            Source list is sorted and collected from start up to PAR.NTASK
        """
        # Apply wildcard rule and check for available sources, exit if no
        # sources found because then we can't proceed
        wildcard = f"{self.source_prefix}_*"
        fids = sorted(glob(os.path.join(PATH.SPECFEM_DATA, wildcard)))
        if not fids:
            print(msg.cli("No matching source files when searching PATH for"
                          "the given WILDCARD",
                          items=[f"PATH: {PATH.SPECFEM_DATA}",
                                 f"WILDCARD: {wildcard}"], header="error"
                          )
                  )
            sys.exit(-1)

        # Create internal definition of sources names by stripping prefixes
        names = [os.path.basename(fid).split("_")[-1] for fid in fids]
        self._source_names = names[:PAR.NTASK]

    def check_solver_parameter_files(self):
        """
        Optional method

        !!! Can be implemented by subclass !!!
        """
        pass

    @property
    def taskid(self):
        """
        Returns the currently running process for embarassingly parallelized
        tasks.

        :rtype: int
        :return: task id for given solver
        """
        return system.taskid()

    @property
    def source_name(self):
        """
        Returns name of source currently under consideration

        :rtype: str
        :return: given source name for given task id
        """
        return self.source_names[self.taskid]

    @property
    def cwd(self):
        """
        Returns working directory currently in use

        :rtype: str
        :return: current solver working directory
        """
        return os.path.join(PATH.SOLVER, self.source_name)

    @property
    def source_names(self):
        """
        Returns list of source names

        :rtype: list
        :return: list of source names
        """
        if self._source_names is None:
            self.check_source_names()

        return self._source_names

    @property
    def mesh_properties(self):
        """
        Returns mesh properties

        :rtype: Struct
        :return: Structure of mesh properties
        """
        if self._mesh_properties is None:
            self.check_mesh_properties()

        return self._mesh_properties

    @property
    def data_filenames(self):
        """
        Template filenames for accessing data

        !!! Must be implemented by subclass !!!
        """
        return NotImplementedError

    @property
    def model_databases(self):
        """
        Template filenames for accessing models

        !!! Must be implemented by subclass !!!
        """
        return NotImplementedError

    @property
    def kernel_databases(self):
        """
        Template filenames for accessing kernels

        !!! Must be implemented by subclass !!!
        """
        return NotImplementedError

    @property
    def source_prefix(self):
        """
        Preferred source prefix

        :rtype: str
        :return: source prefix
        """
        return PAR.SOURCE_PREFIX.upper()


