#!/usr/bin/env python3
"""
This is the subclass seisflows.solver.Specfem3D
This class provides utilities for the Seisflows solver interactions with
Specfem3D Cartesian. It inherits all attributes from seisflows.solver.Base,
and overwrites these functions to provide specified interaction with Specfem3D
"""
import os
import sys
import logging
from glob import glob

import seisflows.plugins.solver.specfem3d as solvertools
from seisflows.tools import unix, msg
from seisflows.tools.wrappers import exists
from seisflows.config import custom_import, SeisFlowsPathsParameters
from seisflows.tools.specfem import getpar, setpar


# Seisflows configuration
PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]

system = sys.modules["seisflows_system"]
preprocess = sys.modules["seisflows_preprocess"]


class Specfem3D(custom_import("solver", "base")):
    """
    Python interface to Specfem3D Cartesian. This subclass inherits functions
    from seisflows.solver.Base

    !!! See base class for method descriptions !!!
    """
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type logger: Logger
        :param logger: Class-specific logging module, log statements pushed
            from this logger will be tagged by its specific module/classname
        """
        super().__init__()

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        # Define the Parameters required by this module
        sf.par("NT", required=True, par_type=float,
               docstr="Number of time steps set in the SPECFEM Par_file")

        sf.par("DT", required=True, par_type=float,
               docstr="Time step or delta set in the SPECFEM Par_file")

        sf.par("FORMAT", required=True, par_type=float,
               docstr="Format of synthetic waveforms used during workflow, "
                      "available options: ['ascii', 'su']")

        sf.par("SOURCE_PREFIX", required=False, default="CMTSOLUTION",
               par_type=str,
               docstr="Prefix of SOURCE files in path SPECFEM_DATA. Available "
                      "['CMTSOLUTION', FORCESOLUTION']")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

        acceptable_formats = ["SU", "ASCII"]
        if PAR.FORMAT.upper() not in acceptable_formats:
            raise Exception(f"'FORMAT' must be {acceptable_formats}")

    def generate_data(self, **model_kwargs):
        """
        Generates data using the True model, exports traces to `traces/obs`

        :param model_kwargs: keyword arguments to pass to `generate_mesh`
        """
        # Create the mesh
        self.generate_mesh(**model_kwargs)

        # Run the Forward simulation
        unix.cd(self.cwd)
        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")
        if PAR.ATTENUATION:
            setpar(key="ATTENUATION", val=".true.", file="DATA/Par_file")
        else:
            setpar(key="ATTENUATION", val=".false.", file="DATA/Par_file")

        self.call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xspecfem3D")

        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard)),
                dst=os.path.join("traces", "obs"))

        # Export traces to disk for permanent storage
        if PAR.SAVETRACES:
            self.export_traces(os.path.join(PATH.OUTPUT, "traces", "obs"))

    def generate_mesh(self, model_path, model_name, model_type=None):
        """
        Performs meshing with internal mesher Meshfem3D and database generation

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

        unix.cd(self.cwd)

        # Run mesh generation
        if model_type == "gll":
            self.check_mesh_properties(model_path)

            src = glob(os.path.join(model_path, "*"))
            dst = self.model_databases
            unix.cp(src, dst)

            self.call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xmeshfem3D")
            self.call_solver(mpiexec=PAR.MPIEXEC, 
                             executable="bin/xgenerate_databases")

        # Export the model for future use in the workflow
        if self.taskid == 0:
            self.export_model(os.path.join(PATH.OUTPUT, model_name))

    def eval_func(self, *args, **kwargs):
        """
        Call eval_func from Base class
        """
        super().eval_func(*args, **kwargs)

        # Work around SPECFEM3D conflicting name conventions of SU data
        self.rename_data()

    def forward(self, path="traces/syn"):
        """
        Calls SPECFEM3D forward solver, exports solver outputs to traces dir

        :type path: str
        :param path: path to export traces to after completion of simulation
        """
        # Set parameters and run forward simulation
        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")
        if PAR.ATTENUATION:
            setpar(key="ATTENUATION", val=".true.", file="DATA/Par_file")
        else:
            setpar(key="ATTENUATION", val=".false`.", file="DATA/Par_file")

        self.call_solver(mpiexec=PAR.MPIEXEC,
                         executable="bin/xgenerate_databases")
        self.call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xspecfem3D")

        # Find and move output traces, by default to synthetic traces dir
        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard)),
                dst=path)

    def adjoint(self):
        """
        Calls SPECFEM3D adjoint solver, creates the `SEM` folder with adjoint
        traces which is required by the adjoint solver
        """
        setpar(key="SIMULATION_TYPE", val="3", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".false.", file="DATA/Par_file")
        setpar(key="ATTENUATION", val=".false.", file="DATA/Par_file")

        unix.rm("SEM")
        unix.ln("traces/adj", "SEM")

        self.call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xspecfem3D")

    def check_solver_parameter_files(self):
        """
        Checks solver parameters 
        """
        # Check the number of steps in the SPECFEM2D Par_file
        nt_str, nt, nt_i = getpar(key="NSTEP", file="DATA/Par_file")
        if int(nt) != PAR.NT:
            if self.taskid == 0:
                print(msg.cli(f"SPECFEM3D {nt_str}=={nt} is not equal "
                              f"SeisFlows PAR.NT=={PAR.NT}. Please ensure "
                              f"that these values match in both files.",
                              header="parameter match error", border="=")
                      )
                sys.exit(-1)

        dt_str, dt, dt_i = getpar(key="DT", file="DATA/Par_file")
        if float(dt) != PAR.DT:
            if self.taskid == 0:
                print(msg.cli(f"SPECFEM3D {dt_str}=={dt} is not equal "
                              f"SeisFlows PAR.DT=={PAR.DT}. Please ensure "
                              f"that these values match in both files.",
                              header="parameter match error", border="=")
                      )
                sys.exit(-1)

        # Ensure that NPROC matches the MESH values
        nproc = self.mesh_properties.nproc
        if nproc != PAR.NPROC:
            if self.taskid == 0:
                print(msg.cli(f"SPECFEM3D mesh NPROC=={nproc} is not equal "
                              f"SeisFlows PAR.NPROC=={PAR.NPROC}. "
                              f"Please check that your mesh matches this val.",
                              header="parameter match error", border="=")
                      )
                sys.exit(-1)

        if "MULTIPLES" in PAR:
            raise NotImplementedError

    def initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM

        .. note::
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        # Initialize adjoint traces as zeroes for all data_filenames
        # write to `traces/adj`
        super().initialize_adjoint_traces()

        # Rename data to work around Specfem naming convetions
        self.rename_data()

        # Workaround for Specfem3D's requirement that all components exist,
        # even ones not in use as adjoint traces
        if PAR.FORMAT.upper() == "SU":
            unix.cd(os.path.join(self.cwd, "traces", "adj"))

            for iproc in range(PAR.NPROC):
                for channel in ["x", "y", "z"]:
                    dst = f"{iproc:d}_d{channel}_SU.adj"
                    if not exists(dst):
                        src = f"{iproc:d}_d{PAR.COMPONENTS[0]}_SU.adj"
                        unix.cp(src, dst)

    def rename_data(self):
        """
        Works around conflicting data filename conventions

        Specfem3D's uses different name conventions for regular traces
        and 'adjoint' traces
        """
        if PAR.FORMAT.upper() == "SU":
            files = glob(os.path.join(self.cwd, "traces", "adj", "*SU"))
            unix.rename(old='_SU', new='_SU.adj', names=files)

    @property
    def data_wildcard(self):
        """
        Returns a wildcard identifier for synthetic data

        :rtype: str
        :return: wildcard identifier for channels
        """
        if PAR.FORMAT.upper() == "SU":
            return f"*_d?_SU"
        elif PAR.FORMAT.upper() == "ASCII":
            return f"*.?X?.sem?"

    @property
    def data_filenames(self):
        """
        Returns the filenames of all data, either by the requested components
        or by all available files in the directory.

        :rtype: list
        :return: list of data filenames
        """
        unix.cd(os.path.join(self.cwd, "traces", "obs"))

        if PAR.COMPONENTS:
            components = PAR.COMPONENTS

            if PAR.FORMAT.upper() == "SU":
                return sorted(glob(f"*_d[{components.lower()}]_SU"))
            elif PAR.FORMAT.upper() == "ASCII":
                return sorted(glob(f"*.?X[{components.upper()}].sem?"))
        else:
            if PAR.FORMAT.upper() == "SU":
                return sorted(glob("*_d?_SU"))
            elif PAR.FORMAT.upper() == "ASCII":
                return sorted(glob("*.???.sem?"))

    @property
    def kernel_databases(self):
        """
        The location of databases for kernel outputs, relative to the current
        working directory. 
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
        return PAR.SOURCE_PREFIX.upper()

