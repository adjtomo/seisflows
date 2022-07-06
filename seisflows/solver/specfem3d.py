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
from seisflows.tools.specfem import setpar, getpar


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
            files = glob(self.data_wildcard(comp=self.par.COMPONENTS.lower()))
        else:
            files = glob(self.data_wildcard(comp="?"))
        return sorted(files)

    @property
    def model_databases(self):
        """
        The location of databases for model outputs, usually
        OUTPUT_FILES/DATABASES_MPI. Value is grabbed from the Par_file
        """
        local_path = getpar(key="LOCAL_PATH",
                            file=os.path.join(self.cwd, "DATA", "Par_file"))[1]
        return os.path.join(self.cwd, local_path)

    @property
    def kernel_databases(self):
        """
        The location of databases for kernel outputs, usually the same as
        'model_databases'
        """
        return self.model_databases

    def eval_func(self, path, write_residuals=True):
        """
        Performs forward simulations and evaluates the misfit function using
        the preprocess module. Overrides to add a data renaming call

        .. note::
            This task should be run in parallel by system.run()

        :type path: str
        :param path: directory from which model is imported and where residuals
            will be exported
        :type write_residuals: bool
        :param write_residuals: calculate and export residuals        """
        super().eval_func(path=path, write_residuals=write_residuals)

        # Work around SPECFEM3D conflicting name conventions of SU data
        self._rename_data()

    def _forward(self, output_path):
        """
        Calls SPECFEM3D forward solver, exports solver outputs to traces dir

        :type output_path: str
        :param output_path: path to export traces to after completion of
            simulation expected values are either 'traces/obs' for 'observation'
            data (i.e., synthetics generated by the TRUE model), or
            'traces/syn', for synthetics generated during function evaluations
        """
        # Set parameters and run forward simulation
        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")
        if self.par.ATTENUATION:
            setpar(key="ATTENUATION", val=".true.", file="DATA/Par_file")
        else:
            setpar(key="ATTENUATION", val=".false`.", file="DATA/Par_file")

        self._call_solver(executable="bin/xgenerate_databases",
                          output="fwd_mesher.log")
        self._call_solver(executable="bin/xmeshfem3D", output="fwd_solver.log")

        # Find and move output traces, by default to synthetic traces dir
        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard())),
                dst=output_path)

    def _adjoint(self):
        """
        Calls SPECFEM3D adjoint solver, creates the `SEM` folder with adjoint
        traces which is required by the adjoint solver
        """
        setpar(key="SIMULATION_TYPE", val="3", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".false.", file="DATA/Par_file")

        # Attenuation should always be OFF during adjoint simulations, else
        # you will get a floating point error
        setpar(key="ATTENUATION", val=".false.", file="DATA/Par_file")

        unix.rm("SEM")
        unix.ln("traces/adj", "SEM")

        self._call_solver(executable="bin/xspecfem3D", output="adj_solver.log")

    def _initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM

        .. note::
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        # Initialize adjoint traces as zeroes for all data_filenames
        # write to `traces/adj`
        super()._initialize_adjoint_traces()

        # Rename data to work around Specfem naming convetions
        self._rename_data()

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


