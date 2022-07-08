#!/usr/bin/env python3
"""
This Solver module is in charge of interacting with external numerical solvers
such as SPECFEM (2D/3D/3D_GLOBE). This SPECFEM base class provides general
functions that work with all versions of SPECFEM. Subclasses will provide
additional capabilities unique to each version of SPECFEM.
"""
import os
import sys
import subprocess
from glob import glob

from seisflows.core import Base
from seisflows.plugins import solver_io
from seisflows.tools import msg, unix
from seisflows.tools.specfem import getpar, Model


class Specfem(Base):
    """
    This base class provides an interface through which solver simulations can
    be set up and run. It should not be used by itself, but rather it is meant
    to provide the foundation for the following child classes:

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

    combine, smooth

        Utilities for combining and smoothing kernels

    generate_data, generate_mesh, eval_fwd, forward, adjoint,
    data_filenames, model_databases, kernel_databases, source_prefix

        !!! Required functions which must be implemented by subclass !!!

    """
    def __init__(self):
        """
        These parameters should not be set by the User_!
        Attributes are just initialized as NoneTypes for clarity and docstrings

        :type parameters: list of str
        :param parameters: a list detailing the parameters to be used to
            define the model, available: ['vp', 'vs', 'rho']
        :type _source_names: hidden attribute,
        :param _source_names: the names of all the sources that are being used
            by the solver
        :type logger: Logger
        :param logger: Class-specific logging module, log statements pushed
            from this logger will be tagged by its specific module/classname
        """
        super().__init__()

        self.required.par(
            "MATERIALS", required=False, par_type=str, default="ELASTIC",
            docstr="Material parameters used to define model. Available: "
                   "['ELASTIC': Vp, Vs, 'ACOUSTIC': Vp, 'ISOTROPIC', "
                   "'ANISOTROPIC']"
        )
        self.required.par(
            "DENSITY", required=False, par_type=str, default="CONSTANT",
            docstr="How to treat density during inversion. Available: "
                   "['CONSTANT': Do not update density, "
                   "'VARIABLE': Update density]"
        )
        self.required.par(
            "ATTENUATION", required=False, par_type=bool, default=False,
            docstr="If True, turn on attenuation during forward "
                   "simulations, otherwise set attenuation off. Attenuation "
                   "is always off for adjoint simulations."
        )
        self.required.par(
            "FORMAT", required=False, par_type=float, default="ASCII",
            docstr="Format of synthetic waveforms used during workflow, "
                   "available options: ['ASCII', 'SU']"
        )
        self.required.par(
            "COMPONENTS", required=False, default="ZNE", par_type=str,
            docstr="Components used to generate data, formatted as a single "
                   "string, e.g. ZNE or NZ or E"
        )
        self.required.par(
            "SOLVERIO", required=False, default="fortran_binary", par_type=str,
            docstr="The format external solver files. Available: "
                   "['fortran_binary']"
        )
        self.required.par(
            "SMOOTH_H", required=False, default=0., par_type=float,
            docstr="Gaussian half-width for horizontal smoothing in units of "
                   "meters. If 0., no smoothing applied"
        )
        self.required.par(
            "SMOOTH_V", required=False, default=0., par_type=float,
            docstr="Gaussian half-width for vertical smoothing in units of "
                   "meters"
        )
        self.required.par(
            "TASKTIME_SMOOTH", required=False, default=1, par_type=int,
            docstr="Large radii smoothing may take longer than normal tasks. "
                   "Allocate additional smoothing task time as a multiple of "
                   "TASKTIME"
        )
        self.required.path(
            "SOLVER", required=False,
            default=os.path.join(self.path.WORKDIR, "scratch", "solver"),
            docstr="scratch path to hold solver working directories"
        )
        self.required.path(
            "DATA", required=False,
            docstr="path to a directory containing any external data "
                   "required by the workflow. Catch all directory that "
                   "can be accessed by all modules"
        )
        self.required.path(
            "SPECFEM_BIN", required=False,
            default=os.path.join(self.path.WORKDIR, "specfem", "bin"),
            docstr="path to the SPECFEM binary executables"
        )
        self.required.path(
            "SPECFEM_DATA", required=False,
            default=os.path.join(self.path.WORKDIR, "specfem", "DATA"),
            docstr="path to the SPECFEM DATA/ directory containing the "
                   "'Par_file', 'STATIONS' file and 'CMTSOLUTION' files"
        )
        # Define the Paths required by this module
        self.required.path(
            "MASK", required=False, docstr="Directory to mask files for "
                                           "gradient masking"
        )

        self.parameters = []
        self._source_names = None

    @property
    def _io(self):
        """
        Solver IO module set by User. Located in seisflows.plugins.solver_io

        :rtype: module
        :return: module containing solver input/output options
        """
        return getattr(solver_io, self.par.SOLVERIO)

    @property
    def taskid(self):
        """
        Returns the currently running process for embarassingly parallelized
        tasks.

        :rtype: int
        :return: task id for given solver
        """
        return self.module("system").taskid()

    @property
    def source_names(self):
        """
        Returns list of source names which should be stored in PAR.SPECFEM_DATA
        Source names are expected to match the following wildcard,
        'PREFIX_*' where PREFIX is something like 'CMTSOLUTION' or 'FORCE'

        :rtype: list
        :return: list of source names
        """
        if self._source_names is None:
            self._check_source_names()
        return self._source_names

    @property
    def source_name(self):
        """
        Returns name of source currently under consideration

        :rtype: str
        :return: given source name for given task id
        """
        return self.source_names[self.taskid]

    @property
    def source_prefix(self):
        """
        Preferred source prefix

        TODO remove this and replace all instances with the self.par parameter

        :rtype: str
        :return: source prefix
        """
        return self.par.SOURCE_PREFIX.upper()

    @property
    def cwd(self):
        """
        Returns working directory currently in use by a running solver instance

        :rtype: str
        :return: current solver working directory
        """
        return os.path.join(self.path.SOLVER, self.source_name)

    def data_wildcard(self, comp="?"):
        """
        Provide a wildcard string that will match the name of the output
        synthetic seismograms

        :type comp: str
        :param comp: single letter defining the component that can be inserted
            into the wildcard. Defaults to '?'
        :rtype: str
        :return: a wildcard string that can be used to search for data
        """
        return NotImplementedError

    @property
    def data_filenames(self):
        """
        A list of waveform files matching the `data_wildcard` which is used
        to keep track of data

        :rtype: list
        :return: list of data filenames
        """
        return NotImplementedError

    @property
    def model_databases(self):
        """
        SPECFEM directory where model database files are saved. This directory
        is SPECFEM version dep.
        """
        return NotImplementedError

    @property
    def kernel_databases(self):
        """
        SPECFEM directory where kernel database files are saved. This directory
        is SPECFEM version dep.
        """
        return NotImplementedError

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

        # Check that other modules have set parameters that will be used here
        for required_parameter in ["NPROC"]:
            assert (required_parameter in self.par), \
                f"Solver requires {required_parameter}"

        available_materials = ["ELASTIC", "ACOUSTIC",  # specfem2d, specfem3d
                               "ISOTROPIC", "ANISOTROPIC"]  # specfem3d_globe
        assert(self.par.MATERIALS.upper() in available_materials), \
            f"MATERIALS must be in {available_materials}"

        acceptable_densities = ["CONSTANT", "VARIABLE"]
        assert(self.par.DENSITY.upper() in acceptable_densities), \
            f"DENSITY must be in {acceptable_densities}"

        acceptable_formats = ["SU", "ASCII"]
        if self.par.FORMAT.upper() not in acceptable_formats:
            raise Exception(f"'FORMAT' must be {acceptable_formats}")

        # Internal parameter list based on user-input material choices
        # Important to reset parameters to a blank list and let the check
        # statements fill it. If not, each time workflow is resumed, parameters
        # list will append redundant parameters and things stop working
        self.parameters = []
        if self.par.MATERIALS.upper() == "ELASTIC":
            assert(self.par.SOLVER.lower() in ["specfem2d", "specfem3d"])
            self.parameters += ["vp", "vs"]
        elif self.par.MATERIALS.upper() == "ACOUSTIC":
            assert(self.par.SOLVER.lower() in ["specfem2d", "specfem3d"])
            self.parameters += ["vp"]
        elif self.par.MATERIALS.upper() == "ISOTROPIC":
            assert(self.par.SOLVER.lower() in ["specfem3d_globe"])
            self.parameters += ["vp", "vs"]
        elif self.par.MATERIALS.upper() == "ANISOTROPIC":
            assert(self.par.SOLVER.lower() in ["specfem3d_globe"])
            self.parameters += ["vpv", "vph", "vsv", "vsh", "eta"]

        if self.par.DENSITY.upper() == "VARIABLE":
            self.parameters.append("rho")

        assert hasattr(solver_io, self.par.SOLVERIO)
        assert hasattr(self._io, "read_slice"), "IO method has no attribute 'read'"
        assert hasattr(self._io, "write_slice"), "IO method has no attribute 'write'"

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
        self._initialize_solver_directories()
        self.generate_data()
        self._set_model(model_name="init", model_type="gll")
        self._initialize_adjoint_traces()

    def _set_model(self, model_name, model_type=None):
        """
        Mesh and database files should have been created during the manual set
        up phase. This function simply checks the mesh properties of that mesh
        and ensures that it is locatable by future SeisFlows processes.

        :type model_name: str
        :param model_name: name of the model to be used as identification
        :type model_type: str
        :param model_type: available model types to be passed to the Specfem3D
            Par_file. See Specfem3D Par_file for available options.
        """
        unix.cd(self.cwd)
        available_model_types = ["gll"]

        # Check the type of model. So far SeisFlows only accepts GLL models
        model_type = model_type or getpar(key="MODEL", file="DATA/Par_file")[1]
        assert(model_type in available_model_types), \
            f"{model_type} not in available types {available_model_types}"

        # Determine which model will be set as the starting model
        if model_name.upper() == "INIT":
            model_path = self.path.MODEL_INIT
        elif model_name.upper() == "TRUE":
            model_path = self.path.MODEL_TRUE
        else:
            raise ValueError(f"model name must be 'INIT' or 'TRUE'")
        assert(os.path.exists(model_path)), f"model {model_path} does not exist"

        if model_type == "gll":
            # Copy the model files (ex: proc000023_vp.bin ...) into database dir
            src = glob(os.path.join(model_path, "*"))
            dst = self.model_databases
            unix.cp(src, dst)

        # Export the model into output folder, ready to be used by other tasks
        if self.taskid == 0:
            self._export_model(os.path.join(self.path.OUTPUT, model_name))

    def generate_data(self):
        """
        Generates observation data to be compared to synthetics. This must
        only be run once. If `PAR.CASE`=='data', then real data will be copied
        over.

        If `PAR.CASE`=='synthetic' then the external solver will use the
        True model to generate 'observed' synthetics. Finally exports traces to
        'cwd/traces/obs'

        Elif `PAR.CASE`=='DATA', will look in PATH.DATA for directories matching
        the given source name and copy ANY files that exist there. e.g., if
        source name is '001', you must store waveform data in PATH.DATA/001/*

        Also exports observed data to OUTPUT if desired
        """
        # If synthetic inversion, generate 'data' with solver
        if self.par.CASE.upper() == "SYNTHETIC":
            if self.path.MODEL_TRUE is not None:
                if self.taskid == 0:
                    self.logger.info("generating 'data' with MODEL_TRUE")
                # Generate synthetic data on the fly using the true model
                self._set_model(model_name="true", model_type="gll")
                self._forward(output_path=os.path.join("traces", "obs"))
        # If Data provided by user, copy directly into the solver directory
        elif self.path.DATA is not None and os.path.exists(self.path.DATA):
            unix.cp(
                src=glob(os.path.join(self.path.DATA, self.source_name, "*")),
                dst=os.path.join("traces", "obs")
            )
        # Save observation data to disk
        if self.par.SAVETRACES:
            self._export_traces(
                path=os.path.join(self.path.OUTPUT, "traces", "obs")
            )

    def eval_func(self, path, write_residuals=True):
        """
        Performs forward simulations and evaluates the misfit function using
        the preprocess module.

        .. note::
            This task should be run in parallel by system.run()

        :type path: str
        :param path: directory from which model is imported and where residuals
            will be exported
        :type write_residuals: bool
        :param write_residuals: calculate and export residuals
        """
        preprocess = self.module("preprocess")

        if self.taskid == 0:
            self.logger.info("running forward simulations")

        unix.cd(self.cwd)
        self._import_model(path)
        self._forward(output_path=os.path.join("traces", "syn"))

        if write_residuals:
            if self.taskid == 0:
                self.logger.debug("calling preprocess.prepare_eval_grad()")
            preprocess.prepare_eval_grad(cwd=self.cwd, taskid=self.taskid,
                                         source_name=self.source_name,
                                         filenames=self.data_filenames
                                         )
            self._export_residuals(path)

    def eval_grad(self, path, export_traces=False):
        """
        Evaluates gradient by carrying out adjoint simulations.

        .. note::
            It is expected that eval_func() has already been run as this
            function looks for adjoint sources in 'cwd/traces/adj'

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

        self._adjoint()
        self._export_kernels(path)

        if export_traces:
            self._export_traces(path=os.path.join(path, "traces", "syn"),
                                prefix="traces/syn")
            self._export_traces(path=os.path.join(path, "traces", "adj"),
                                prefix="traces/adj")

    def postprocess_kernels(self, path_grad):
        """
        Sums kernels from individual sources, with optional smoothing

        .. note::
            This function needs to be run on system, i.e., called by
            system.run(single=True)

        :type path_grad: str
        :param path_grad: directory containing sensitivity kernels in the
            scratch directory to be summed and smoothed. Output summed and
            summed + smoothed kernels will be saved here as well.
        """
        # If specified, smooth the kernels in the vertical and horizontal and
        # save both (summed, summed+smoothed) to separate output directories
        kernel_path = os.path.join(path_grad, "kernels")
        path_sum_nosmooth = os.path.join(kernel_path, "sum_nosmooth")
        path_sum = os.path.join(kernel_path, "sum")

        if (self.par.SMOOTH_H > 0) or (self.par.SMOOTH_V > 0):
            self.logger.debug(f"saving un-smoothed and summed kernels to:\n"
                              f"{path_sum_nosmooth}")
            self.combine(input_path=kernel_path, output_path=path_sum_nosmooth)

            self.logger.info(f"smoothing gradient: H={self.par.SMOOTH_H}m, "
                             f"V={self.par.SMOOTH_V}m")
            self.logger.debug(f"saving smoothed kernels to:\n{path_sum}")
            self.smooth(input_path=path_sum_nosmooth, output_path=path_sum,
                          span_h=self.par.SMOOTH_H, span_v=self.par.SMOOTH_V)

        # Combine all the input kernels, generating the unscaled gradient
        else:
            self.logger.debug(f"saving summed kernels to:\n{path_sum}")
            self.combine(input_path=kernel_path, output_path=path_sum)


    # def apply_hess(self, path):
    #     """
    #     High level solver interface that computes action of Hessian on a given
    #     model vector. A gradient evaluation must have already been carried out.
    #
    #     TODO preprocess has no function prepare_apply_hess()
    #
    #     :type path: str
    #     :param path: directory to which output files are exported
    #     """
    #     raise NotImplementedError("must be implemented by solver subclass")
    #
    #     unix.cd(self.cwd)
    #     self.import_model(path)
    #     unix.mkdir("traces/lcg")
    #
    #     self.forward("traces/lcg")
    #     preprocess.prepare_apply_hess(self.cwd)
    #     self.adjoint()
    #     self.export_kernels(path)

    def _forward(self, output_path):
        """
        Calls forward solver with the appropriate parameters in the Par_file set
        Also exports data to the correct output_path

        :type output_path: str
        :param output_path: path to export traces to after completion of
            simulation expected values are either 'traces/obs' for 'observation'
            data (i.e., synthetics generated by the TRUE model), or
            'traces/syn', for synthetics generated during function evaluations
        """
        raise NotImplementedError("must be implemented by solver subclass")

    def _adjoint(self):
        """
        Calls adjoint solver with the appropriate parameters in the Par_file set
        Also takes care of setting up the SEM/ directory where SPECFEM expects
        adjoint sources to be
        """
        raise NotImplementedError("must be implemented by solver subclass")

    def _call_solver(self, executable, output="solver.log"):
        """
        Calls MPI solver executable to run solver binaries, used by individual
        processes to run the solver on system. If the external solver returns a
        non-zero exit code (failure), this function will return a negative
        boolean.

        :type executable: str
        :param executable: executable function to call
        :type output: str
        :param output: where to redirect stdout
        """
        # Executable may come with additional sub arguments, we only need to
        # check that the actually executable exists
        if not os.path.exists(executable.split(" ")[0]):
            print(msg.cli(f"solver executable {executable} does not exist",
                          header="external solver error", border="="))
            sys.exit(-1)

        # mpiexec is None when running in serial mode, so e.g., ./xmeshfem2D
        if self.par.SYSTEM in ["workstation"]:
            exc_cmd = f"./{executable}"
        # Otherwise mpiexec is system dependent (e.g., srun, mpirun)
        else:
            exc_cmd = f"{self.par.MPIEXEC} {executable}"

        # Run solver. Write solver stdout (log files) to text file
        try:
            with open(output, "w") as f:
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
        unix.cd(self.cwd)

        if parameters is None:
            parameters = self.parameters

        if not os.path.exists(output_path):
            unix.mkdir(output_path)

        # Write the source names into the kernel paths file for SEM/ directory
        with open("kernel_paths", "w") as f:
            f.writelines(
                [os.path.join(input_path, f"{name}\n")
                 for name in self.source_names]
            )

        # Call on xcombine_sem to combine kernels into a single file
        for name in parameters:
            # e.g.: mpiexec ./bin/xcombine_sem alpha_kernel kernel_paths output
            self._call_solver(executable=" ".join([f"bin/xcombine_sem",
                                                   f"{name}_kernel",
                                                   "kernel_paths",
                                                   output_path]
                                                  )
                              )

    def smooth(self, input_path, output_path, parameters=None, span_h=0.,
               span_v=0., output="smooth.log"):
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
        unix.cd(self.cwd)

        if parameters is None:
            parameters = self.parameters

        if not os.path.exists(output_path):
            unix.mkdir(output_path)

        # mpiexec ./bin/xsmooth_sem SMOOTH_H SMOOTH_V name input output use_gpu
        for name in parameters:
            self._call_solver(executable=" ".join(["bin/xsmooth_sem",
                                                   str(span_h), str(span_v),
                                                   f"{name}_kernel",
                                                   os.path.join(input_path, ""),
                                                   os.path.join(output_path, ""),
                                                   ".false"]),
                              output=output)

        # Rename output files
        files = glob(os.path.join(output_path, "*"))
        unix.rename(old="_smooth", new="", names=files)

    def _import_model(self, path):
        """
        File transfer utility to move a SPEFEM2D model into the correct location
        for a workflow.

        :type path: str
        :param path: path to the SPECFEM2D model
        :return:
        """
        unix.cp(src=glob(os.path.join(path, "model", "*")),
                dst=os.path.join(self.cwd, "DATA")
                )

    def _import_traces(self, path):
        """
        File transfer utility. Import traces into the workflow.

        :type path: str
        :param path: path to traces
        """
        src = glob(os.path.join(path, 'traces', self.source_name, '*'))
        dst = os.path.join(self.cwd, 'traces', 'obs')
        unix.cp(src, dst)

    def _export_model(self, path, parameters=None):
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

    def _export_kernels(self, path):
        """
        File transfer utility. Export kernels to disk

        :type path: str
        :param path: path to save kernels
        """
        if self.taskid == 0:
            self.logger.debug(f"exporting kernels to:\n{path}")

        unix.cd(self.kernel_databases)

        # Work around conflicting name conventions
        self._rename_kernels()

        src = glob("*_kernel.bin")
        dst = os.path.join(path, "kernels", self.source_name)
        unix.mkdir(dst)
        unix.mv(src, dst)

    def _export_residuals(self, path):
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

    def _export_traces(self, path, prefix="traces/obs"):
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

    def _rename_kernels(self):
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

    def _initialize_solver_directories(self):
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
            self.logger.info(f"initializing {self.par.NTASK} solver directories")

        unix.rm(self.cwd)
        unix.mkdir(self.cwd)
        unix.cd(self.cwd)

        # Create directory structure
        for cwd_dir in ["bin", "DATA", "OUTPUT_FILES/DATABASES_MPI",
                        "traces/obs", "traces/syn", "traces/adj",
                        self.model_databases, self.kernel_databases]:
            unix.mkdir(cwd_dir)

        # Copy exectuables into the bin/ directory
        src = glob(os.path.join(self.path.SPECFEM_BIN, "*"))
        dst = os.path.join("bin", "")
        unix.cp(src, dst)

        # Copy all input files except source files
        src = glob(os.path.join(self.path.SPECFEM_DATA, "*"))
        src = [_ for _ in src if self.source_prefix not in _]
        dst = os.path.join("DATA", "")
        unix.cp(src, dst)

        # Symlink event source specifically, strip the source name as SPECFEM
        # just expects `source_name`
        src = os.path.join(self.path.SPECFEM_DATA,
                           f"{self.source_prefix}_{self.source_name}")
        dst = os.path.join("DATA", self.source_prefix)
        unix.ln(src, dst)

        if self.taskid == 0:
            mainsolver = os.path.join(self.path.SOLVER, "mainsolver")
            # Symlink taskid_0 as mainsolver in solver directory for convenience
            if not os.path.exists(mainsolver):
                unix.ln(self.cwd, mainsolver)
                self.logger.debug(f"source {self.source_name} symlinked as "
                                  f"mainsolver")
        else:
            # Copy the initial model from mainsolver into current directory
            # Avoids the need to run multiple instances of xgenerate_databases
            src = os.path.join(self.path.SOLVER, "mainsolver", "OUTPUT_FILES",
                               "DATABASES_MPI")
            dst = os.path.join(self.cwd, "OUTPUT_FILES", "DATABASES_MPI")
            unix.cp(src, dst)

    def _initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM.
        This is only done for the 'base' the Preprocess class.

        .. note::
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        preprocess = self.module("preprocess")

        if self.par.PREPROCESS.upper() == "DEFAULT":
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

    def _check_source_names(self):
        """
        Determines names of sources by applying wildcard rule to
        user-supplied input files

        .. note::
            Source list is sorted and collected from start up to PAR.NTASK
        """
        # Apply wildcard rule and check for available sources, exit if no
        # sources found because then we can't proceed
        wildcard = f"{self.source_prefix}_*"
        fids = sorted(glob(os.path.join(self.path.SPECFEM_DATA, wildcard)))
        if not fids:
            print(msg.cli("No matching source files when searching PATH for"
                          "the given WILDCARD",
                          items=[f"PATH: {self.path.SPECFEM_DATA}",
                                 f"WILDCARD: {wildcard}"], header="error"
                          )
                  )
            sys.exit(-1)

        # Create internal definition of sources names by stripping prefixes
        names = [os.path.basename(fid).split("_")[-1] for fid in fids]
        self._source_names = names[:self.par.NTASK]



