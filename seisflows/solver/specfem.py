#!/usr/bin/env python3
"""
This Solver module is in charge of interacting with external numerical solvers
such as SPECFEM (2D/3D/3D_GLOBE). This SPECFEM base class provides general
functions that work with all versions of SPECFEM. Subclasses will provide
additional capabilities unique to each version of SPECFEM.

.. note::
    The Base class implementation is almost completely SPECFEM2D related.
    However, SPECFEM2D requires a few unique parameters that 3D/3D_GLOBE
    do not. Because of the inheritance architecture of SeisFlows, we do not
    want the 3D and 3D_GLOBE versions to inherit 2D-specific parameters, so
    we need this this more generalized SPECFEM base class.

TODO
    - add in `apply_hess` functionality that was partially written in legacy code
    - move `_initialize_adjoint_traces` to workflow.migration
    - Add density scaling based on Vp?
"""
import os
import sys
import subprocess
from glob import glob

from seisflows import logger
from seisflows.tools import msg, unix
from seisflows.tools.config import get_task_id, Dict
from seisflows.tools.specfem import getpar, setpar, check_source_names


class Specfem:
    """
    Solver SPECFEM
    --------------
    Defines foundational structure for Specfem-based solver module. 
    Generalized SPECFEM interface to manipulate SPECFEM2D/3D/3D_GLOBE w/ Python

    Parameters
    ----------
    :type syn_data_format: str
    :param syn_data_format: data format for reading synthetic traces into memory.
        Available: ['SU': seismic unix format, 'ASCII': human-readable ascii]
    :type materials: str
    :param materials: Material parameters used to define model. Available:
        ['ELASTIC': Vp, Vs, 'ACOUSTIC': Vp, 'ISOTROPIC', 'ANISOTROPIC']
    :type density: bool
    :param density: How to treat density during inversion. If True, updates
        density during inversion. If False, keeps it constant.
        TODO allow density scaling during an inversion
    :type attenuation: bool
    :param attenuation: How to treat attenuation during inversion.
        if True, turns on attenuation during forward simulations only. If
        False, attenuation is always set to False. Requires underlying
        attenution (Q_mu, Q_kappa) model
    :type smooth_h: float
    :param smooth_h: Gaussian half-width for horizontal smoothing in units
        of meters. If 0., no smoothing applied. Only applicable for workflows:
        ['migration', 'inversion'], ignored for 'forward' workflow.
        SPECFEM3D_GLOBE only: if `smooth_type`=='laplacian' then this is just 
        the X and Y extent of the applied smoothing
    :type smooth_h: float
    :param smooth_v: Gaussian half-width for vertical smoothing in units
        of meters. Only applicable for workflows: ['migration', 'inversion'],
        ignored for 'forward' workflow.
        SPECFEM3D_GLOBE only: if `smooth_type`=='laplacian' then this is just 
        the Z extent of the applied smoothing
    :type components: str
    :param components: components to consider and tag data with. Should be
        string of letters such as 'RTZ'
    :type solver_io: str
    :param solver_io: format of model/kernel/gradient files expected by the
        numerical solver. Available: ['fortran_binary': default .bin files].
        TODO: ['adios': ADIOS formatted files]
    :type source_prefix: str
    :param source_prefix: prefix of source/event/earthquake files. If None,
        will attempt to guess based on the specific solver chosen.
    :type mpiexec: str
    :param mpiexec: MPI executable used to run parallel processes. Should also
        be defined for the system module

    Paths
    -----
    :type path_data: str
    :param path_data: path to any externally stored data required by the solver
    :type path_specfem_bin: str
    :param path_specfem_bin: path to SPECFEM bin/ directory which
        contains binary executables for running SPECFEM
    :type path_specfem_data: str
    :param path_specfem_data: path to SPECFEM DATA/ directory which must
        contain the CMTSOLUTION, STATIONS and Par_file files used for
        running SPECFEM
    ***
    """
    def __init__(self, syn_data_format="ascii",  materials="acoustic",
                 density=False, nproc=1, ntask=1, attenuation=False,
                 smooth_h=0., smooth_v=0., components="ZNE",
                 source_prefix=None, mpiexec=None, workdir=os.getcwd(),
                 path_solver=None, path_eval_grad=None,
                 path_data=None, path_specfem_bin=None, path_specfem_data=None,
                 path_model_init=None, path_model_true=None, path_output=None,
                 **kwargs):
        """
        Set default SPECFEM interface parameters

        .. note::
            Paths listed here are shared with `workflow.forward` and so are not
            included in the class docstring.

        :type workdir: str
        :param workdir: working directory in which to look for data and store
            results. Defaults to current working directory
        :type path_solver: str
        :param path_solver: scratch path for all solver related tasks
        :type path_model_init: str
        :param path_model_init: path to the starting model used to calculate the
            initial misfit. Must match the expected `solver_io` format.
        :type path_model_true: str
        :param path_model_true: path to a target model if `case`=='synthetic' and
            a set of synthetic 'observations' are required for workflow.
        :type path_output: str
        :param path_output: shared output directory on disk for more permanent
            storage of solver related files such as traces, kernels, gradients.
        """
        # Publically accessible parameters
        self.syn_data_format = syn_data_format
        self.materials = materials
        self.nproc = nproc
        self.ntask = ntask
        self.density = density
        self.attenuation = attenuation
        self.smooth_h = smooth_h
        self.smooth_v = smooth_v
        self.components = components
        self.source_prefix = source_prefix or "SOURCE"

        # Define internally used directory structure
        self.path = Dict(
            scratch=path_solver or os.path.join(workdir, "scratch", "solver"),
            eval_grad=path_eval_grad or
                      os.path.join(workdir, "scratch", "evalgrad"),
            data=path_data or os.path.join(workdir, "SFDATA"),
            output=path_output or os.path.join(workdir, "output"),
            specfem_bin=path_specfem_bin,
            specfem_data=path_specfem_data,
            model_init=path_model_init,
            model_true=path_model_true,
        )
        self.path.mainsolver = os.path.join(self.path.scratch, "mainsolver")

        # Private internal parameters for keeping track of solver requirements
        self._parameters = []
        if self.density:
            self._parameters.append("rho")

        self._mpiexec = mpiexec
        self._source_names = None  # for property source_names
        self._ext = ""  # for database file extensions

        # Define available choices for check parameters
        self._available_model_types = ["gll"]
        self._available_materials = [
            "ELASTIC", "ACOUSTIC",  # specfem2d, specfem3d
            "ISOTROPIC", "ANISOTROPIC"  # specfem3d_globe
        ]
        # SPECFEM2D specific attributes. Should be overwritten by 3D versions
        self._syn_available_data_formats = ["ASCII", "SU"]
        self._required_binaries = ["xspecfem2D", "xmeshfem2D", "xcombine_sem",
                                   "xsmooth_sem"]
        self._acceptable_source_prefixes = ["SOURCE", "FORCE", "FORCESOLUTION"]
        
        # Empty variable that will need to be overwritten by SPECFEM3D_GLOBE
        self._regions = None

    def check(self):
        """
        Checks parameter validity for SPECFEM input files and model parameters
        """
        assert(self.materials.upper() in self._available_materials), \
            f"solver.materials must be in {self._available_materials}"

        if self.syn_data_format.upper() not in self._syn_available_data_formats:
            raise NotImplementedError(
                f"solver.syn_data_format must be {self._syn_available_data_formats}"
            )

        # Check that User has provided appropriate binary files to run SPECFEM
        assert(self.path.specfem_bin is not None and
               os.path.exists(self.path.specfem_bin)), (
            f"`path_specfem_bin` must exist and must point to directory " 
            f"containing SPECFEM executables"
        )
        for fid in self._required_binaries:
            assert(os.path.exists(os.path.join(self.path.specfem_bin, fid))), (
                f"`path_specfem_bin`/{fid} does not exist but is required by "
                f"SeisFlows solver module"
            )

        # Check that SPECFEM/DATA directory exists
        assert(self.path.specfem_data is not None and
               os.path.exists(self.path.specfem_data)), (
            f"`path_specfem_data` must exist and must point to directory " 
            f"containing SPECFEM input files"
        )
        for fid in ["STATIONS", "Par_file"]:
            assert(os.path.exists(os.path.join(self.path.specfem_data, fid))), (
                f"DATA/{fid} does not exist but is required by SeisFlows solver"
            )

        # Make sure source files exist and are appropriately labeled
        assert(self.source_prefix in self._acceptable_source_prefixes), (
            f"SPECFEM `source_prefix` must be in "
            f"{self._acceptable_source_prefixes}"
            )
        assert(glob(os.path.join(self.path.specfem_data,
                                 f"{self.source_prefix}*"))), (
            f"No source files with prefix {self.source_prefix} found in DATA/")

        # Check that model type is set correctly in the Par_file
        model_type = getpar(key="MODEL",
                            file=os.path.join(self.path.specfem_data,
                                              "Par_file"))[1]
        assert(model_type in self._available_model_types), (
            f"SPECFEM Par_file parameter `model`='{model_type}' does not "
            f"match acceptable model types: {self._available_model_types}"
            )

        # Assign file extensions to be used for database file searching
        if model_type == "gll":
            self._ext = ".bin"

        # Make sure the initial model is set and actually contains files
        assert(self.path.model_init is not None and
               os.path.exists(self.path.model_init)), \
            f"`path_model_init` is required for the solver, but does not exist"

        assert(len(glob(os.path.join(self.path.model_init, "*")))), \
            f"`path_model_init` is empty but should have model files"

        if self.path.model_true is not None:
            assert(os.path.exists(self.path.model_true)), \
                f"`path_model_true` is provided but does not exist"
            assert(len(glob(os.path.join(self.path.model_true, "*")))), \
                f"`path_model_true` is empty but should have model files"

        # Check that the number of tasks/events matches the number of events
        self._source_names = check_source_names(
            path_specfem_data=self.path.specfem_data,
            source_prefix=self.source_prefix, ntask=self.ntask
        )

        assert(isinstance(self.density, bool)), \
            f"solver `density` must be True (variable) or False (constant)"

        # Check the size of the DATA/ directory and let the User know if 
        # large files are present, e.g., tomo xyz files or topo/bathy
        for root, dirs, files in os.walk(self.path.specfem_data):
            for name in files:
                fullpath = os.path.join(root, name)
                if not os.path.islink(fullpath):
                    filesize = os.path.getsize(fullpath) / 1E9  # Bytes -> GB
                    if filesize > 0.5:
                        logger.warning(
                            f"SPECFEM DATA/ file '{fullpath}' is >.5GB and "
                            f"will be copied {self.ntask} time(s). Please be "
                            f"sure to check if this file is necessary for your "
                            f"workflow"
                            )

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
            self._source_names = check_source_names(
                path_specfem_data=self.path.specfem_data,
                source_prefix=self.source_prefix, ntask=self.ntask
            )
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
        return self.source_names[get_task_id()]

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
        if self.syn_data_format.upper() == "SU":
            return f"U{comp}_file_single.su"
        elif self.syn_data_format.upper() == "ASCII":
            return f"*.?X{comp}.sem?"

    def model_wildcard(self, par="*", kernel=False):
        """
        Returns a wildcard identifier to search for models kernels generated by
        the solver. An example SPECFEM2D/3D kernel filename (in 
        FORTRAN binary file format) is: 'proc000001_rho_kernel.bin'
        Whereas the corresponding model would be 'proc000001_rho.bin'

        Allows dynamically searching for specific files when renaming, moving
        or copying files. Also allows for different wildcard for 3D_GLOBE 
        version

        :type comp: str
        :param comp: component formatter, defaults to wildcard '?'
        :type kernel: bool
        :param kernel: wildcarding a kernel file. If True, adds the 'kernel' 
            tag. If not, assuming we are wildcarding for a model file
        :rtype: str
        :return: wildcard identifier for channels
        """
        if kernel:
            _ker = "_kernel"
        else:
            _ker = ""
        return f"proc??????_{par}{_ker}{self._ext}"

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

        if self.components:
            comp_glob = f"[{self.components}]"  # e.g., [NEZ]
        else:
            comp_glob = "?"
        data_wildcard = self.data_wildcard(comp=comp_glob)
        file_glob = os.path.join(self.cwd, "traces", choice, data_wildcard)
        filenames = glob(file_glob)

        if not filenames:
            logger.critical(
                msg.cli("The property `solver.data_filenames`, used to search "
                        "for waveform files, is empty and should not be. "
                        "Please check solver parameters: ",
                        items=[f"failed wildcard: {file_glob}"],
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
        :return: path where SPECFEM2D database files are stored, relative to
            `solver.cwd`
        """
        return "DATA"

    @property
    def model_files(self):
        """
        Return a list of paths to model files that match the internal parameter
        list. Used to generate model vectors of the same length as gradients.

        :rtype: list
        :return: a list of full paths to model files that matches the internal
            list of solver parameters
        """
        _model_files = []
        for par in self._parameters:
            _model_files += glob(os.path.join(self.path.mainsolver,
                                              self.model_databases,
                                              self.model_wildcard(par=par))
                                              )
        return _model_files

    @property
    def kernel_databases(self):
        """
        The location of kernel inputs and outputs as defined by SPECFEM2D
        This is RELATIVE to a SPECFEM2D working directory.

         .. note::
            This path is SPECFEM version dependent so SPECFEM3D/3D_GLOBE
            versions must overwrite this function

        :rtype: str
        :return: path where SPECFEM2D database files are stored, relative to
            `solver.cwd`
        """
        return "OUTPUT_FILES"

    def setup(self):
        """
        Prepares solver scratch directories for an impending workflow.

        Sets up directory structure expected by SPECFEM and copies or generates
        seismic data to be inverted or migrated.

        Exports INIT/STARTING and TRUE/TARGET models to disk (output/ dir.)
        """
        self._initialize_working_directories()
        self._export_starting_models()

    def forward_simulation(self, executables=None, save_traces=False,
                           export_traces=False, **kwargs):
        """
        Wrapper for SPECFEM binaries: 'xmeshfem?D' 'xgenerate_databases',
                                      'xspecfem?D'

        Calls SPECFEM2D forward solver, exports solver outputs to traces dir

         .. note::
            SPECFEM3D/3D_GLOBE versions must overwrite this function

        :type executables: list or None
        :param executables: list of SPECFEM executables to run, in order, to
            complete a forward simulation. This can be left None in most cases,
            which will select default values based on the specific solver
            being called (2D/3D/3D_GLOBE). It is made an optional parameter
            to keep the function more general for inheritance purposes.
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
        if executables is None:
            executables = ["bin/xmeshfem2D", "bin/xspecfem2D"]

        unix.cd(self.cwd)
        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")

        # Calling subprocess.run() for each of the binary executables listed
        for exc in executables:
            # e.g., fwd_mesher.log
            stdout = f"fwd_{self._exc2log(exc)}.log"
            self._run_binary(executable=exc, stdout=stdout)

        # Work around SPECFEM's version dependent file names
        if self.syn_data_format.upper() == "SU":
            for tag in ["d", "v", "a", "p"]:
                unix.rename(old=f"single_{tag}.su", new="single.su",
                            names=glob(os.path.join("OUTPUT_FILES", "*.su")))
        # Exporting traces to disk (output/) for more permanent storage
        if export_traces:
            if not os.path.exists(export_traces):
                unix.mkdir(export_traces)
            unix.cp(
                src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard())),
                dst=export_traces
            )
        # Save traces somewhere else in the scratch/ directory for easier access
        if save_traces:
            if not os.path.exists(save_traces):
                unix.mkdir(save_traces)
            unix.mv(
                src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard())),
                dst=save_traces
            )

    def adjoint_simulation(self, executables=None, save_kernels=False,
                           export_kernels=False):
        """
        Wrapper for SPECFEM binary 'xspecfem?D'

        Calls SPECFEM2D adjoint solver, creates the `SEM` folder with adjoint
        traces which is required by the adjoint solver. Renames kernels
        after they have been created from 'alpha' and 'beta' to 'vp' and 'vs',
        respectively.

         .. note::
            SPECFEM3D/3D_GLOBE versions must overwrite this function

        :type executables: list or None
        :param executables: list of SPECFEM executables to run, in order, to
            complete an adjoint simulation. This can be left None in most cases,
            which will select default values based on the specific solver
            being called (2D/3D/3D_GLOBE). It is made an optional parameter
            to keep the function more general for inheritance purposes.
        :type save_kernels: str
        :param save_kernels: move the kernels from their native SPECFEM output
            location to another path. This is used to move kernels to another
            SeisFlows scratch directory so that they are discoverable by
            other modules. The typical location they are moved to is
            path_eval_grad
        :type export_kernels: str
        :param export_kernels: export/copy/save kernels from the scratch
            directory to a more permanent storage location. i.e., copy files
            from their original location. Note that kernel file sizes are LARGE,
            so exporting kernels can lead to massive storage requirements.
        """
        if executables is None:
            executables = ["bin/xspecfem2D"]

        unix.cd(self.cwd)

        setpar(key="SIMULATION_TYPE", val="3", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".false.", file="DATA/Par_file")

        unix.rm("SEM")
        unix.ln("traces/adj", "SEM")

        # Calling subprocess.run() for each of the binary executables listed
        for exc in executables:
            # e.g., adj_solver.log
            stdout = f"adj_{self._exc2log(exc)}.log"
            logger.info(f"running SPECFEM executable {exc}, log to '{stdout}'")
            self._run_binary(executable=exc, stdout=stdout)

        # Rename kernels to work w/ conflicting name conventions
        # Change directory so that the rename doesn't affect the full path
        # Deals with both SPECFEM3D and 3D_GLOBE, which adds in the 'reg?' tag
        unix.cd(self.kernel_databases)
        for tag in ["alpha", "alpha[hv]", "reg?_alpha", "reg?_alpha[hv]"]:
            names = glob(self.model_wildcard(par=tag, kernel=True))
            if names:
                logger.debug(f"renaming output event kernels: '{tag}' -> 'vp'")
                unix.rename(old="alpha", new="vp", names=names)

        for tag in ["beta", "beta[hv]", "reg?_beta", "reg?_beta[hv]"]:
            names = glob(self.model_wildcard(par=tag, kernel=True))
            if names:
                logger.debug(f"renaming output event kernels: '{tag}' -> 'vs'")
                unix.rename(old="beta", new="vs", names=names)

        # Save and export the kernels to user-defined locations
        if export_kernels:
            unix.mkdir(export_kernels)
            for par in self._parameters:
                unix.cp(src=glob(self.model_wildcard(par=par, kernel=True)),
                        dst=export_kernels)

        if save_kernels:
            unix.mkdir(save_kernels)
            for par in self._parameters:
                unix.mv(src=glob(self.model_wildcard(par=par, kernel=True)),
                        dst=save_kernels)

    def combine(self, input_path, output_path, parameters=None):
        """
        Wrapper for 'xcombine_sem'.
        Sums kernels from individual source contributions to create gradient.

        .. note::
            The binary xcombine_sem simply sums matching databases

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
            defaults to `self._parameters`
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
            # e.g., smooth_vp.log
            stdout = f"{self._exc2log(exc)}_{name}.log"
            self._run_binary(executable=exc, stdout=stdout)

    def smooth(self, input_path, output_path, parameters=None, span_h=None,
               span_v=None, use_gpu=False):
        """
        Wrapper for SPECFEM binary: xsmooth_sem
        Smooths kernels by convolving them with a 3D Gaussian

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
            defaults to `self._parameters`
        :type span_h: float
        :param span_h: horizontal smoothing length in meters
        :type span_v: float
        :param span_v: vertical smoothing length in meters
        :type use_gpu: bool
        :param use_gpu: whether to use GPU acceleration for smoothing. Requires
            GPU compiled binaries and GPU compute node.
        """
        unix.cd(self.cwd)

        # Assign some default parameters from class attributes if not given
        if parameters is None:
            parameters = self._parameters
        if span_h is None:
            span_h = self.smooth_h
        if span_v is None:
            span_v = self.smooth_v

        logger.debug(f"smoothing {parameters} with horizontal Gaussian "
                     f"{span_h}m and vertical Gaussian {span_v}m")

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
            # e.g., combine_vs.log
            stdout = f"{self._exc2log(exc)}_{name}.log"
            self._run_binary(executable=exc, stdout=stdout)

        # Rename output files to remove the '_smooth' suffix which SeisFlows
        # will not recognize
        files = glob(os.path.join(output_path, "*"))
        unix.rename(old="_smooth", new="", names=files)

    def _run_binary(self, executable, stdout="solver.log", with_mpi=True):
        """
        Calls MPI solver executable to run solver binaries, used by individual
        processes to run the solver on system. If the external solver returns a
        non-zero exit code (failure), this function will return a negative
        boolean.

        .. note::
            This function ASSUMES it is being run from a SPECFEM working
            directory, i.e., that the executables are located in ./bin/

        .. note::
            This is essentially an error-catching wrapper of subprocess.run()

        :type executable: str
        :param executable: executable function to call. May or may not start
            E.g., acceptable calls for the solver would './bin/xspecfem2D'.
            Also accepts additional command line arguments such as:
            'xcombine_sem alpha_kernel kernel_paths...'
        :type stdout: str
        :param stdout: where to redirect stdout
        :type with_mpi: bool
        :param with_mpi: If `mpiexec` is given, use MPI to run the executable.
            Some executables (e.g., combine_vol_data_vtk) must be run in
            serial so this flag allows them to turn off MPI running.
        :raises SystemExit: If external numerical solver return any failure
            code while running
        """
        # Executable may come with additional sub arguments, we only need to
        # check that the actually executable exists
        if not unix.which(executable.split(" ")[0]):
            logger.critical(msg.cli(f"executable '{executable}' does not exist",
                            header="external solver error", border="="))
            sys.exit(-1)

        # Prepend with `mpiexec` if we are running with MPI
        # looks something like: `mpirun -n 4 ./bin/xspecfem2d`
        if self._mpiexec and with_mpi:
            executable = f"{self._mpiexec} -n {self.nproc} {executable}"
        logger.debug(f"running executable with cmd: '{executable}'")

        try:
            with open(stdout, "w") as f:
                subprocess.run(executable, shell=True, check=True, stdout=f,
                               stderr=f)
        except (subprocess.CalledProcessError, OSError) as e:
            logger.critical(
                msg.cli("The external numerical solver has returned a "
                        "nonzero exit code (failure). Consider stopping any "
                        "currently running jobs to avoid wasted "
                        "computational resources. Check 'scratch/solver/"
                        f"mainsolver/{stdout}' for the solvers stdout log "
                        "message. The failing command and error message are:",
                        items=[f"exc: {executable}", f"err: {e}"],
                        header="external solver error",
                        border="=")
            )
            sys.exit(-1)

    @staticmethod
    def _exc2log(exc):
        """
        Very simple conversion utility to get log file names based on binaries.
        e.g., binary 'xspecfem2D' will return 'solver'. Helps keep log file
        naming consistent and generalizable

        TODO add a check here to see if the log file exists, and then use
            `number_fid` to increment so that we keep all the output logs

        :type exc: str
        :param exc: specfem executable, e.g., xspecfem2D, xgenerate_databases
        :rtype: str
        :return: logfile name that matches executable name
        """
        convert_dict = {"specfem": "solver", "meshfem": "mesher",
                        "generate_databases": "mesher", "smooth": "smooth",
                        "combine": "combine"}
        for key, val in convert_dict.items():
            if key in exc:
                return val
        else:
            return "logger"

    def import_model(self, path_model):
        """
        Copy files from given `path_model` into the current working directory
        model database. Used for grabbing starting models (e.g., MODEL_INIT)
        and models that have been perturbed by the optimization library.

        :type path_model: str
        :param path_model: path to an existing starting model
        """
        assert(os.path.exists(path_model)), f"model {path_model} does not exist"
        unix.cd(self.cwd)

        # Copy the model files (ex: proc000023_vp.bin ...) into database dir
        src = glob(os.path.join(path_model, f"*{self._ext}"))
        dst = os.path.join(self.cwd, self.model_databases, "")
        unix.cp(src, dst)

    def _initialize_working_directories(self):
        """
        Serial task used to initialize working directories for each of the a
        available sources
        """
        logger.info(f"initializing {self.ntask} solver directories")
        for source_name in self.source_names:
            cwd = os.path.join(self.path.scratch, source_name)
            if os.path.exists(cwd):
                continue
            self._initialize_working_directory(cwd=cwd)

    def _initialize_working_directory(self, cwd=None):
        """
        Creates scratch directory structure expected by SPECFEM
        (i.e., bin, DATA, OUTPUT_FILES). Copies executables (bin) and
        input data (DATA) directories, prepares simulation input files.

        Each directory will act as completely independent Specfem working dir.
        This allows for embarrassing parallelization while avoiding the need
        for intra-directory communications, at the cost of temporary disk space.

        .. note::
            path to binary executables must be supplied by user as SeisFlows has
            no mechanism for automatically compiling from source code.

        :type cwd: str
        :param cwd: optional scratch working directory to intialize. If None,
            will set based on current running seisflows task (self.taskid)
        """
        # Define a constant list of required SPECFEM dir structure, relative cwd
        _required_structure = {"bin", "DATA", "OUTPUT_FILES", "traces/obs", 
                               "traces/syn", "traces/adj", self.model_databases,
                               self.kernel_databases}

        # Allow this function to be called on system or in serial
        if cwd is None:
            cwd = self.cwd
            source_name = self.source_name
        else:
            source_name = os.path.basename(cwd)

        logger.debug(f"initializing solver directory source: {source_name}")

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
        dst = os.path.join(cwd, "DATA", "")
        unix.cp(src, dst)

        # Symlink event source specifically, only retain source prefix
        src = os.path.join(self.path.specfem_data,
                           f"{self.source_prefix}_{source_name}")
        dst = os.path.join(cwd, "DATA", self.source_prefix)
        unix.ln(src, dst)

        # Symlink TaskID==0 as mainsolver in solver directory for convenience
        if self.source_names.index(source_name) == 0:
            if not os.path.exists(self.path.mainsolver):
                logger.debug(f"linking source '{source_name}' as 'mainsolver'")
                unix.ln(cwd, self.path.mainsolver)

    def _export_starting_models(self, parameters=None):
        """
        Export the initial and target models to the SeisFlows output/ directory.

        :type parameters: list
        :param parameters: list of parameters to export. If None, will default
            to `self._parameters`
        """
        if parameters is None:
            parameters = self._parameters

        # Export the initial and target models to the SeisFlows output directory
        for name, model in zip(["MODEL_INIT", "MODEL_TRUE"],
                               [self.path.model_init, self.path.model_true]):
            # Skip over if user has not provided model path (e.g., real data
            # inversion will not have `model_true`)
            if not model:
                continue
            dst = os.path.join(self.path.output, name, "")
            if not os.path.exists(dst):
                unix.mkdir(dst)
            for par in parameters:
                src = glob(os.path.join(model, f"*{par}{self._ext}"))
                unix.cp(src, dst)
