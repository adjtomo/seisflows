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

    @property
    def _io(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self._io

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

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

    def setup(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self.setup()

    def _set_model(self, model_name, model_type="gll"):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._set_model(model_name=model_name, model_type=model_type)

    def generate_data(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self.generate_data()

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

    def eval_grad(self, path, export_traces=False):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self.eval_grad(path=path, export_traces=export_traces)

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
        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard)),
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

    def _call_solver(self, executable, output="solver.log"):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._call_solver(executable=executable, output=output)

    def load(self, path, prefix="", suffix="", parameters=None):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.load(path=path, prefix=prefix, suffix=suffix,
                         parameters=parameters)

    def save(self, save_dict,  path, parameters=None, prefix="", suffix=""):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self.save(save_dict=save_dict, path=path, parameters=parameters,
                  prefix=prefix, suffix=suffix)

    def merge(self, model, parameters=None):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.merge(model=model, parameters=parameters)

    def split(self, m, parameters=None):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.split(m=m, parameters=parameters)

    def combine(self, input_path, output_path, parameters=None):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.combine(input_path=input_path, output_path=output_path,
                            parameters=parameters)

    def smooth(self, input_path, output_path, parameters=None, span_h=0.,
               span_v=0., output="smooth.log"):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self.smooth(input_path=input_path, output_path=output_path,
                           parameters=parameters, span_h=span_h,
                           span_v=span_v, output=output)

    def _import_model(self, path):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self._import_model(path=path)

    def _import_traces(self, path):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self._import_traces(path=path)

    def _export_model(self, path, parameters=None):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self._export_model(path=path, parameters=parameters)

    def _export_kernels(self, path):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self._export_kernels(path=path)

    def _export_residuals(self, path):
        """Inherits from seisflows.solver.specfem.Specfem"""
        return self._export_residuals(path=path)

    def _export_traces(self, path, prefix="traces/obs"):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._export_traces(path=path, prefix=prefix)

    def _rename_kernels(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._rename_kernels()

    def _initialize_solver_directories(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._initialize_solver_directories()

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

    def _check_mesh_properties(self, path=None):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._check_mesh_properties(path=path)

    def _check_source_names(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._check_source_names()

    def _rename_data(self):
        """
        Works around conflicting data filename conventions

        Specfem3D's uses different name conventions for regular traces
        and 'adjoint' traces
        """
        if self.par.FORMAT.upper() == "SU":
            files = glob(os.path.join(self.cwd, "traces", "adj", "*SU"))
            unix.rename(old='_SU', new='_SU.adj', names=files)


