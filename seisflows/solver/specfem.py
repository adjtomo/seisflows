#!/usr/bin/env python3
"""
This Solver module is in charge of interacting with external numerical solvers
such as SPECFEM (2D/3D/3D_GLOBE). This SPECFEM base class provides general
functions that work with all versions of SPECFEM. Subclasses will provide
additional capabilities unique to each version of SPECFEM.

TODO add in `apply_hess` functionality that was partially written in legacy code
TODO move `_initialize_adjoint_traces` to workflow.migration
"""
import os
import sys
import subprocess
from glob import glob

from seisflows import logger
from seisflows.core import Dict
from seisflows.plugins import solver_io as solver_io_dir
from seisflows.tools import msg, unix
from seisflows.tools.specfem import getpar, setpar


class Specfem:
    """
    This base class provides an interface through which solver simulations can
    be set up and run. It should not be used by itself, but rather it is meant
    to provide the foundation for: SPECFEM2D/3D/3D_GLOBE

    .. note::
        This Base class implementation is almost completely SPECFEM2D related.
        However, SPECFEM2D requires a few unique parameters that 3D/3D_GLOBE
        do not. Because of the inheritance architecture of SeisFlows, we do not
        want the 3D and 3D_GLOBE versions to inherit 2D-specific parameters, so
        we need this this more generalized SPECFEM base class.

    .. note:::
        This class supports only acoustic and isotropic elastic inversions.
    """
    def __init__(self, case="data", data_format="ascii",  materials="elastic",
                 density=False, nproc=1, ntask=1, attenuation=False,
                 components="ZNE", solver_io="fortran_binary",
                 source_prefix=None,mpiexec=None, path_solver=None,
                 path_data=None, path_specfem_bin=None,
                 path_specfem_data=None, path_model_init=None,
                 path_model_true=None, path_output=None, **kwargs):
        """
        SPECFEM Solver parameters

        :type data_format: str
        :param data_format: data format for reading traces into memory.
            Availalble: ['SU' seismic unix format, 'ASCII' human-readable ascii]
        :type materials: str
        :param materials: Material parameters used to define model. Available:
            ['ELASTIC': Vp, Vs, 'ACOUSTIC': Vp, 'ISOTROPIC', 'ANISOTROPIC']
        :type density: bool
        :param density: How to treat density during inversion. If True, updates
            density during inversion. If False, keeps it constant.
            TODO Add density scaling based on Vp?
        :type attenuation: bool
        :param attenuation: How to treat attenuation during inversion.
            if True, turns on attenuation during forward simulations only. If
            False, attenuation is always set to False. Requires underlying
            attenution (Q_mu, Q_kappa) model
        :type components: str
        :param components: components to consider and tag data with. Should be
            string of letters such as 'RTZ'
        :type path.specfem_bin: str
        :param path.specfem_bin: path to SPECFEM bin/ directory which
            contains binary executables for running SPECFEM
        :type path.specfem_data: str
        :param path.specfem_data: path to SPECFEM DATA/ directory which must
            contain the CMTSOLUTION, STATIONS and Par_file files used for
            running SPECFEM

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
        self.case = case
        self.data_format = data_format
        self.materials = materials
        self.nproc = nproc
        self.ntask = ntask
        self.density = density
        self.attenuation = attenuation
        self.components = components
        self.solver_io = solver_io
        self.mpiexec = mpiexec

        # Define internally used directory structure
        _cwd = os.getcwd()
        self.path = Dict(
            scratch=path_solver or os.path.join(_cwd, "scratch", "solver"),
            data=path_data or os.path.join(_cwd, "SFDATA"),
            output=path_output or os.path.join(_cwd, "output"),
            mainsolver=os.path.join(_cwd, "scratch", "mainsolver"),
            specfem_bin=path_specfem_bin,
            specfem_data=path_specfem_data,
            model_init=path_model_init,
            model_true=path_model_true,
        )

        # Establish internally defined parameter system
        self._parameters = []
        if self.density:
            self._parameters.append("rho")

        self._source_names = None
        self._io = getattr(solver_io_dir, self.solver_io)

        # Define available choices for check parameter
        self._available_model_types = ["gll"]
        self._available_materials = [
            "ELASTIC", "ACOUSTIC",  # specfem2d, specfem3d
            "ISOTROPIC", "ANISOTROPIC"  # specfem3d_globe
        ]
        self._available_data_formats = ["ASCII", "SU"]
        self._required_binaries = ["xspecfem2d", "xmeshfem2d", "xcombine_sem",
                                   "xsmooth_sem"]

        # These are parameters that need to be established by child classes
        self.source_prefix = source_prefix
        self._acceptable_source_prefixes = []

    @property
    def taskid(self):
        """
        Returns the currently running process for embarassingly parallelized
        tasks. Task IDs are assigned to the environment by system.run().
        Task IDs are simply integer values from 0 to the number of
        simultaneously running tasks.

        .. note::
            Dependent on environment variable 'SEISFLOWS_TASKID' which is
            assigned by system.run() to each individually running process.

        :rtype: int
        :return: task id for given solver
        """
        _taskid = os.getenv("SEISFLOWS_TASKID")
        if _taskid is None:
            _taskid = 0
            logger.warning("Environment variable 'SEISFLOWS_TASKID' not found. "
                           "Assigning Task ID == 0")
        return int(_taskid)

    @property
    def source_names(self):
        """
        Returns list of source names which should be stored in PAR.SPECFEM_DATA
        Source names are expected to match the following wildcard,
        'PREFIX_*' where PREFIX is something like 'CMTSOLUTION' or 'FORCE'

        .. note::
            Dependent on environment variable 'SEISFLOWS_TASKID' which is
            assigned by system.run() to each individually running process.

        :rtype: list
        :return: list of source names
        """
        if self._source_names is None:
            self._source_names = self._check_source_names()
        return self._source_names

    @property
    def source_name(self):
        """
        Returns name of source currently under consideration

        .. note::
            Dependent on environment variable 'SEISFLOWS_TASKID' which is
            assigned by system.run() to each individually running process.

        :rtype: str
        :return: given source name for given task id
        """
        return self.source_names[self.taskid]

    @property
    def cwd(self):
        """
        Returns working directory currently in use by a running solver instance

        .. note::
            Dependent on environment variable 'SEISFLOWS_TASKID' which is
            assigned by system.run() to each individually running process.

        :rtype: str
        :return: current solver working directory
        """
        return os.path.join(self.path.scratch, self.source_name)

    def data_wildcard(self, comp="?"):
        """
        Returns a wildcard identifier for synthetic data based on SPECFEM2D
        file naming schema. Allows formatting dcomponent e.g.,
        when called by solver.data_filenames.

        .. note::
            SPECFEM3D/3D_GLOBE versions must overwrite this function

        :type comp: str
        :param comp: component formatter, defaults to wildcard '?'
        :rtype: str
        :return: wildcard identifier for channels
        """
        if self.data_format.upper() == "SU":
            return f"U{comp}_file_single.su"
        elif self.data_format.upper() == "ASCII":
            return f"*.?X{comp}.sem?"

    def data_filenames(self, choice="obs"):
        """
        Returns the filenames of SPECFEM2D data, either by the requested
        components or by all available files in the directory.

         .. note::
            SPECFEM3D/3D_GLOBE versions must overwrite this function

        .. note::
            If the glob returns an  empty list, this function exits the
            workflow because filenames should not be empty is they're being
            queried

        :rtype: list
        :return: list of data filenames
        """

        assert(choice in ["obs", "syn", "adj"]), \
            f"choice must be: 'obs', 'syn' or 'adj'"
        unix.cd(os.path.join(self.cwd, "traces", choice))

        filenames = []
        if self.components:
            for comp in self.components:
                filenames = glob(self.data_wildcard(comp=comp.lower()))
        else:
            filenames = glob(self.data_wildcard(comp="?"))

        if not filenames:
            print(msg.cli("The property solver.data_filenames, used to search "
                          "for traces in 'scratch/solver/*/traces' is empty "
                          "and should not be. Please check solver parameters: ",
                          items=[f"data_wildcard: {self.data_wildcard()}"],
                          header="data filenames error", border="=")
                  )
            sys.exit(-1)

        return filenames

    @property
    def model_databases(self):
        """
        The location of model inputs and outputs as defined by SPECFEM2D.
        This is RELATIVE to a SPECFEM2D working directory.

         .. note::
            This path is SPECFEM version dependent so SPECFEM3D/3D_GLOBE
            versions must overwrite this function

        :rtype: str
        :return: path where SPECFEM2D database files are stored
        """
        return "DATA"

    @property
    def kernel_databases(self):
        """
        The location of kernel inputs and outputs as defined by SPECFEM2D
        This is RELATIVE to a SPECFEM2D working directory.

         .. note::
            This path is SPECFEM version dependent so SPECFEM3D/3D_GLOBE
            versions must overwrite this function

        :rtype: str
        :return: path where SPECFEM2D database files are stored
        """
        return "OUTPUT_FILES"

    def check(self):
        """
        Checks parameter validity for SPECFEM input files and model parameters
        """
        assert(self.case.upper() in ["DATA", "SYNTHETIC"]), \
            f"solver.case must be 'DATA' or 'SYNTHETIC'"

        assert(self.materials.upper() in self._available_materials), \
            f"solver.materials must be in {self._available_materials}"

        if self.data_format.upper() not in self._available_data_formats:
            raise NotImplementedError(
                f"solver.data_format must be {self._available_data_formats}"
            )

        # Make sure we can read in the model/kernel/gradient files
        assert hasattr(solver_io_dir, self.solver_io)
        assert hasattr(self._io, "read_slice"), \
            "IO method has no attribute 'read'"
        assert hasattr(self._io, "write_slice"), \
            "IO method has no attribute 'write'"

        # Check that User has provided appropriate bin/ and DATA/ directories
        for name, dir_ in zip(["bin/", "DATA/"],
                              [self.path.specfem_bin, self.path.specfem_data]):
            assert(dir_ is not None), f"SPECFEM path '{name}' cannot be None"
            assert(os.path.exists(dir_)), f"SPECFEM path '{name}' must exist"

        # Check that the required SPECFEM files are available
        for fid in [self.source_prefix, "STATIONS", "Par_file"]:
            assert(os.path.exists(os.path.join(self.path.specfem_data, fid))), (
                f"DATA/{fid} does not exist but is required by SeisFlows solver"
            )

        # Check that required binary files exist which are called upon by solver
        for fid in self._required_binaries:
            assert(os.path.exists(os.path.join(self.path.specfem_bin, fid))), (
                f"bin/{fid} does not exist but is required by SeisFlows solver"
            )

        # Check that the 'case' variable matches required models
        if self.case.upper() == "SYNTHETIC":
            assert(os.path.exists(self.path.model_true)), (
                f"solver.case == 'synthetic' requires `path.model_true`"
            )

        # Make sure source files exist and are appropriately labeled
        assert(self.source_prefix in self._acceptable_source_prefixes)
        assert(glob(os.path.join(self.path.specfem_data,
                                 f"{self.source_prefix}*"))), (
            f"No source files with prefix {self.source_prefix} found in DATA/")

        # Check that model type is set correctly in the Par_file
        model_type = getpar(key="MODEL",
                            file=os.path.join(self.path.specfem_data,
                                              "Par_file"))[1]
        assert(model_type in self._available_model_types), \
            f"{model_type} not in available types {self._available_model_types}"

    def setup(self):
        """
        Prepares solver scratch directories for an impending workflow.

        Sets up directory structure expected by SPECFEM and copies or generates
        seismic data to be inverted or migrated

        TODO the .bin during model export assumes GLL file format, more general?
        """
        self._initialize_working_directories()

        # Export the initial model to the SeisFlows output directory
        unix.mkdir(self.path.output)
        for key in self._parameters:
            src = glob(os.path.join(self.path.model_init, f"*{key}.bin"))
            dst = os.path.join(self.path.output, "MODEL_INIT", "")
            unix.cp(src, dst)

        # TODO move this into workflow.migration
        # self._initialize_adjoint_traces()

    # def generate_data(self, save_traces=False):
    #     """
    #     Generates observation data to be compared to synthetics. This must
    #     only be run once. If `PAR.CASE`=='data', then real data will be copied
    #     over.
    #  TODO move this to workflow
    #
    #     If `PAR.CASE`=='synthetic' then the external solver will use the
    #     True model to generate 'observed' synthetics. Finally exports traces to
    #     'cwd/traces/obs'
    #
    #     Elif `PAR.CASE`=='DATA', will look in PATH.DATA for directories matching
    #     the given source name and copy ANY files that exist there. e.g., if
    #     source name is '001', you must store waveform data in PATH.DATA/001/*
    #
    #     Also exports observed data to OUTPUT if desired
    #     """
    #     # If synthetic inversion, generate 'data' with solver
    #     if self.case.upper() == "SYNTHETIC":
    #         if self.path.model_true is not None:
    #             if self.taskid == 0:
    #                 logger.info("generating 'data' with MODEL_TRUE")
    #
    #             # Generate synthetic data on the fly using the true model
    #             self.import_model(path_model=self.path.model_true)
    #             self.forward_simulation(
    #                 save_traces=os.path.join("traces", "obs")
    #             )
    #     # If Data provided by user, copy directly into the solver directory
    #     elif self.path.data is not None and os.path.exists(self.path.data):
    #         unix.cp(
    #             src=glob(os.path.join(self.path.data, self.source_name, "*")),
    #             dst=os.path.join(self.cwd, "traces", "obs")
    #         )
    #
    #     # Save observation data to disk
    #     if save_traces:
    #         self._export_traces(
    #             path=os.path.join(self.path.output, "traces", "obs")
    #         )

    def generate_data(self, export_traces=False):
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

        :type export_traces: str
        :param export_traces: path to copy and save traces to a more permament
            storage location as waveform stored in scratch/ are liable to be
            deleted or overwritten
        """
        # Basic checks to make sure there are True model files to copy
        assert(self.case.upper() == "SYNTHETIC")
        assert(os.path.exists(self.path.model_true))
        assert(glob(os.path.join(self.path.model_true, "*")))

        # Generate synthetic data on the fly using the true model
        self.import_model(path_model=self.path.model_true)
        self.forward_simulation(
            save_traces=os.path.join(self.cwd, "traces", "obs"),
            export_traces=export_traces
        )

    def import_data(self):
        """
        Import data from an existing directory into the current working
        directory, required if 'observed' waveform data will be provided by
        the User rather than automatically collected (with Pyatoa) or generated
        synthetically (with external solver)
        """
        # Simple checks to make sure we can actually import data
        assert(self.case.upper() == "DATA")
        assert(self.path.data is not None)
        assert(os.path.exists(os.path.join(self.path.data, self.source_name)))
        assert(glob(os.path.join(self.path.data, self.source_name, "*")))

        src = os.path.join(self.path.data, self.source_name, "*")
        dst = os.path.join(self.cwd, "traces", "obs")

        unix.cp(src, dst)

    def eval_func(self, path, preprocess=None):
        """
        Performs forward simulations and evaluates the misfit function using
        the preprocess module. solver.eval_func is bundled with
        preprocess.prepare_eval_grad because they are meant to be run serially
        so it is better to lump them together into a single allocation.

        .. note::
            This task should be run in parallel by system.run()

        :type path: str
        :param path: directory from which model is imported and where residuals
            will be exported
        :type preprocess: instance
        :param preprocess: SeisFlows preprocess module which can be used to
            prepare gradient evaluation by comparing misfit and creating
            adjoint sources. If None, only forward simulations will be
            performed
        """
        unix.cd(self.cwd)

        if self.taskid == 0:
            logger.info("running forward simulations")

        self._import_model(path)
        self._forward(output_path=os.path.join("traces", "syn"))

        if preprocess:
            if self.taskid == 0:
                logger.debug("call preprocess to prepare gradient evaluation")
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
            logger.debug("running adjoint simulations")

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

    def forward_simulation(self, save_traces=False, export_traces=False):
        """
        Calls SPECFEM2D forward solver, exports solver outputs to traces dir

         .. note::
            SPECFEM3D/3D_GLOBE versions must overwrite this function

        :type save_traces: str
        :param save_traces: move files from their native SPECFEM output location
            to another directory. This is used to move output waveforms to
            'traces/obs' or 'traces/syn' so that SeisFlows knows where to look
            for them, and so that SPECFEM doesn't overwrite existing files
            during subsequent forward simulations
        :type export_traces: str
        :param export_traces: export traces from the scratch directory to a more
            permanent storage location. i.e., copy files from their original
            location
        """
        unix.cd(self.cwd)

        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")

        self._call_solver(executable="bin/xmeshfem2D", output="fwd_mesher.log")
        self._call_solver(executable="bin/xspecfem2D", output="fwd_solver.log")

        # Work around SPECFEM2D's version dependent file names
        if self.data_format.upper() == "SU":
            for tag in ["d", "v", "a", "p"]:
                unix.rename(old=f"single_{tag}.su", new="single.su",
                            names=glob(os.path.join("OUTPUT_FILES", "*.su")))

        if export_traces:
            unix.cp(
                src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard())),
                dst=export_traces
            )

        if save_traces:
            unix.mv(
                src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard())),
                dst=save_traces
            )

    def adjoint_simulation(self):
        """
        Calls SPECFEM2D adjoint solver, creates the `SEM` folder with adjoint
        traces which is required by the adjoint solver.

         .. note::
            SPECFEM3D/3D_GLOBE versions must overwrite this function
        """
        unix.cd(self.cwd)

        setpar(key="SIMULATION_TYPE", val="3", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".false.", file="DATA/Par_file")

        unix.rm("SEM")
        unix.ln("traces/adj", "SEM")

        # Deal with different SPECFEM2D name conventions for regular traces and
        # "adjoint" traces
        if self.data_format.upper == "SU":
            unix.rename(old=".su", new=".su.adj",
                        names=glob(os.path.join("traces", "adj", "*.su")))

        self._call_solver(executable="bin/xspecfem2D", output="adj_solver.log")

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
        if not self.mpiexec:
            exc_cmd = f"./{executable}"
        # Otherwise mpiexec is system dependent (e.g., srun, mpirun)
        else:
            exc_cmd = f"{self.mpiexec} {executable}"

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
        :type output_path: strs
        :param output_path: path to export the outputs of xcombine_sem
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        """
        unix.cd(self.cwd)

        if parameters is None:
            parameters = self._parameters

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
            # e.g.: mpiexec bin/xcombine_sem alpha_kernel kernel_paths output/
            exc = f"bin/xcombine_sem {name}_kernel kernel_paths {output_path}"
            self._call_solver(executable=exc)

    def smooth(self, input_path, output_path, parameters=None, span_h=0.,
               span_v=0., use_gpu=False, output="smooth.log"):
        """
        Postprocessing wrapper: xsmooth_sem
        Smooths kernels by convolving them with a Gaussian.

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
        :type use_gpu: bool
        :param use_gpu: whether to use GPU acceleration for smoothing. Requires
            GPU compiled binaries and GPU compute node.
        """
        unix.cd(self.cwd)

        if parameters is None:
            parameters = self._parameters

        if not os.path.exists(output_path):
            unix.mkdir(output_path)

        # Ensure trailing '/' character, required by xsmooth_sem
        input_path = os.path.join(input_path, "")
        output_path = os.path.join(output_path, "")
        if use_gpu:
            use_gpu = ".true"
        else:
            use_gpu = ".false"
        # mpiexec ./bin/xsmooth_sem SMOOTH_H SMOOTH_V name input output use_gpu
        for name in parameters:
            exc = (f"bin/xsmooth_sem {str(span_h)} {str(span_v)} {name}_kernel "
                   f"{input_path} {output_path} {use_gpu}")
            self._call_solver(executable=exc, output=output)

        # Rename output files to remove the '_smooth' suffix which SeisFlows
        # will not recognize
        files = glob(os.path.join(output_path, "*"))
        unix.rename(old="_smooth", new="", names=files)

    def import_model(self, path_model):
        """
        Copy files from given `path_model` into the current working directory
        model database

        :type path_model: str
        :param path_model: path to an existing starting model
        """
        assert(os.path.exists(path_model)), f"model {path_model} does not exist"
        unix.cd(self.cwd)

        # Copy the model files (ex: proc000023_vp.bin ...) into database dir
        src = glob(os.path.join(path_model, "*"))
        dst = os.path.join(self.cwd, self.model_databases, "")
        unix.cp(src, dst)

    def _export_model(self):
        """
        File transfer utility. Export the model to disk. Run from master solver.
        """
        unix.mkdir(self.path.output)
        for key in self._parameters:
            files = glob(os.path.join(self.model_databases, f"*{key}.bin"))
            unix.cp(files, self.path.output)

    def _export_kernels(self):
        """
        File transfer utility. Export kernels to disk
        """
        unix.cd(self.kernel_databases)

        # Work around conflicting name conventions
        self._rename_kernels()

        src = glob("*_kernel.bin")
        dst = os.path.join(self.path.output, "kernels", self.source_name)
        unix.mkdir(dst)
        unix.mv(src, dst)

    def _export_residuals(self, path):
        """
        File transfer utility. Export residuals to disk.

        :type path: str
        :param path: path to save residuals
        """
        if self.taskid == 0:
            logger.debug(f"exporting residuals to:\n{path}")

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
            logger.debug("exporting traces to {path} {prefix}")

        unix.mkdir(os.path.join(path))

        src = os.path.join(self.cwd, prefix)
        dst = os.path.join(path, self.source_name)
        unix.cp(src, dst)

    @staticmethod
    def _rename_kernels():
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

    def _initialize_working_directories(self):
        """
        Serial task used to initialize working directories for each of the a
        available sources

        TODO run this with concurrent futures for speedup?
        """
        logger.info(f"initializing {self.ntask} solver directories")
        for source_name in self.source_names:
            cwd = os.path.join(self.path, source_name)
            self._initialize_working_directory(cwd=cwd)

    def _initialize_working_directory(self, cwd=None):
        """
        Creates directory structure expected by SPECFEM
        (i.e., bin/, DATA/, OUTPUT_FILES/). Copies executables and prepares
        input files.

        Each directory will act as completely independent Specfem working dir.
        This allows for embarrassing parallelization while avoiding the need
        for intra-directory communications, at the cost of temporary disk space.

        .. note::
            Path to binary executables must be supplied by user as SeisFlows has
            no mechanism for automatically compiling from source code.

        :type cwd: str
        :param cwd: optional scratch working directory to intialize. If None,
            will set based on current running seisflows task (self.taskid)
        """
        # Define a constant list of required SPECFEM dir structure, relative cwd
        _required_structure = ["bin", "DATA",
                               "traces/obs", "traces/syn", "traces/adj",
                               self.model_databases, self.kernel_databases]

        # Allow this function to be called on system or in serial
        if cwd is None:
            cwd = self.cwd
            _source_name = os.path.basename(cwd)
            taskid = self.source_names.index(_source_name)
        else:
            cwd = self.cwd
            taskid = self.taskid

        if taskid == 0:
            logger.info(f"initializing {self.ntask} solver directories")

        # Starting from a fresh working directory
        unix.rm(cwd)
        unix.mkdir(cwd)
        for dir_ in _required_structure:
            unix.mkdir(os.path.join(cwd, dir_))

        # Copy existing SPECFEM exectuables into the bin/ directory
        src = glob(os.path.join(self.path.specfem_bin, "*"))
        dst = os.path.join(cwd, "bin", "")
        unix.cp(src, dst)

        # Copy all input DATA/ files except the source files
        src = glob(os.path.join(self.path.specfem_data, "*"))
        src = [_ for _ in src if self.source_prefix not in _]
        dst = os.path.join("DATA", "")
        unix.cp(src, dst)

        # Symlink event source specifically, only retain source prefix
        src = os.path.join(self.path.specfem_data,
                           f"{self.source_prefix}_{self.source_name}")
        dst = os.path.join("DATA", self.source_prefix)
        unix.ln(src, dst)

        # Symlink TaskID==0 as mainsolver in solver directory for convenience
        if taskid == 0:
            if not os.path.exists(self.path.mainsolver):
                logger.debug(f"symlink {self.source_name} as 'mainsolver'")
                unix.ln(cwd, self.path.mainsolver)

    def _initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM.
        This is only done for the 'base' the Preprocess class.

        TODO move this into workflow setup

        .. note::
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        preprocess = self.module("preprocess")

        if self.par.PREPROCESS.upper() == "DEFAULT":
            if self.taskid == 0:
                logger.debug(f"intializing {len(self.data_filenames)} "
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
        Determines names of sources by applying wildcard rule to user-supplied
        input files. Source names are only provided up to PAR.NTASK and are
        returned in alphabetical order.

        :rtype: list
        :return: alphabetically ordered list of source names up to PAR.NTASK
        """
        assert(self.path.specfem_data is not None), \
            f"solver source names requires 'solver.path.specefm_data' to exist"
        assert(os.path.exists(self.path.specfem_data)), \
            f"solver source names requires 'solver.path.specfem_data' to exist"

        # Apply wildcard rule and check for available sources, exit if no
        # sources found because then we can't proceed
        wildcard = f"{self.source_prefix}_*"
        fids = sorted(glob(os.path.join(self.path.specfem_data, wildcard)))
        if not fids:
            print(msg.cli("No matching source files when searching PATH for "
                          "the given WILDCARD",
                          items=[f"PATH: {self.path.specfem_data}",
                                 f"WILDCARD: {wildcard}"], header="error"
                          )
                  )
            sys.exit(-1)

        # Create internal definition of sources names by stripping prefixes
        names = [os.path.basename(fid).split("_")[-1] for fid in fids]

        return names[:self.ntask]
