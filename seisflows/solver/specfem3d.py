#!/usr/bin/env python3
"""
This is the subclass seisflows.solver.Specfem3D
This class provides utilities for the Seisflows solver interactions with
Specfem3D Cartesian. It inherits all attributes from seisflows.solver.Base,
and overwrites these functions to provide specified interaction with Specfem3D
"""
import os
from glob import glob

from seisflows.solver.specfem import Specfem
from seisflows.tools import unix
from seisflows.tools.wrappers import exists
from seisflows.tools.specfem import getpar, setpar


class Specfem3D(Specfem):
    """
    Python interface to Specfem3D Cartesian.
    """
    def __init__(self):
        """
        Initiate parameters required for Specfem3D Cartesian
        """
        super().__init__()

        self.required.par(
            "SOURCE_PREFIX", required=False, default="CMTSOLUTION",
            par_type=str,
            docstr="Prefix of SOURCE files in path SPECFEM_DATA. Available "
                   "['CMTSOLUTION', FORCESOLUTION']")

    @property
    def io(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.io

    @property
    def taskid(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.taskid

    @property
    def source_names(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.source_names

    @property
    def source_name(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.source_name

    @property
    def source_prefix(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.source_prefix

    @property
    def cwd(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.cwd

    @property
    def mesh_properties(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.mesh_properties

    @property
    def data_wildcard(self, comp="?"):
        """
        Returns a wildcard identifier for synthetic data

        :rtype: str
        :return: wildcard identifier for channels
        """
        if self.par.FORMAT.upper() == "SU":
            return f"*_d?_SU"
        elif self.par.FORMAT.upper() == "ASCII":
            return f"*.?X{comp}.sem?"

    @property
    def data_filenames(self):
        """
        Returns the filenames of all data, either by the requested components
        or by all available files in the directory.

        :rtype: list
        :return: list of data filenames
        """
        unix.cd(os.path.join(self.cwd, "traces", "obs"))

        if self.par.COMPONENTS:
            components = self.par.COMPONENTS

            if self.par.FORMAT.upper() == "SU":
                return sorted(glob(f"*_d[{components.lower()}]_SU"))
            elif self.par.FORMAT.upper() == "ASCII":
                return sorted(glob(f"*.?X[{components.upper()}].sem?"))
        else:
            if self.par.FORMAT.upper() == "SU":
                return sorted(glob("*_d?_SU"))
            elif self.par.FORMAT.upper() == "ASCII":
                return sorted(glob("*.???.sem?"))

    @property
    def model_databases(self):
        """
        The location of databases for model outputs
        """
        return os.path.join(self.cwd, "OUTPUT_FILES", "DATABASES_MPI")

    @property
    def kernel_databases(self):
        """
        The location of databases for kernel outputs
        """
        return self.model_databases

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

    def set_model(self, model_name, model_type="gll"):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self.set_model(model_name=model_name, model_type=model_type)

    def generate_data(self, model_name, model_type="gll"):
        """
        Generates data using the True model, exports traces to `traces/obs`

        :param model_kwargs: keyword arguments to pass to `generate_mesh`
        """
        # Create the mesh

        # Run the Forward simulation
        unix.cd(self.cwd)
        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")
        if self.par.ATTENUATION:
            setpar(key="ATTENUATION", val=".true.", file="DATA/Par_file")
        else:
            setpar(key="ATTENUATION", val=".false.", file="DATA/Par_file")

        self.call_solver(executable="bin/xspecfem3D", output="true_solver.log")

        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard)),
                dst=os.path.join("traces", "obs"))

        # Export traces to disk for permanent storage
        if self.par.SAVETRACES:
            self.export_traces(os.path.join(self.path.OUTPUT, "traces", "obs"))

    def eval_func(self, *args, **kwargs):
        """
        Call eval_func from Base class
        """
        super().eval_func(*args, **kwargs)

        # Work around SPECFEM3D conflicting name conventions of SU data
        self._rename_data()

    def forward(self, path="traces/syn"):
        """
        Calls SPECFEM3D forward solver, exports solver outputs to traces dir

        :type path: str
        :param path: path to export traces to after completion of simulation
        """
        # Set parameters and run forward simulation
        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")
        if self.par.ATTENUATION:
            setpar(key="ATTENUATION", val=".true.", file="DATA/Par_file")
        else:
            setpar(key="ATTENUATION", val=".false`.", file="DATA/Par_file")

        self.call_solver(executable="bin/xgenerate_databases",
                         output="fwd_mesher.log")
        self.call_solver(executable="bin/xmeshfem3D", output="fwd_solver.log")

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

        self.call_solver(executable="bin/xmeshfem3D", output="adj_solver.log")

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
        if self.par.FORMAT.upper() == "SU":
            unix.cd(os.path.join(self.cwd, "traces", "adj"))

            for iproc in range(self.par.NPROC):
                for channel in ["x", "y", "z"]:
                    dst = f"{iproc:d}_d{channel}_SU.adj"
                    if not exists(dst):
                        src = f"{iproc:d}_d{self.par.COMPONENTS[0]}_SU.adj"
                        unix.cp(src, dst)

    def _rename_data(self):
        """
        Works around conflicting data filename conventions

        Specfem3D's uses different name conventions for regular traces
        and 'adjoint' traces
        """
        if self.par.FORMAT.upper() == "SU":
            files = glob(os.path.join(self.cwd, "traces", "adj", "*SU"))
            unix.rename(old='_SU', new='_SU.adj', names=files)


