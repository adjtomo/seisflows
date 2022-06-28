#!/usr/bin/env python3
"""
This is the subclass seisflows.solver.specfem2d

This class provides utilities for the Seisflows solver interactions with
Specfem2D. It inherits all attributes from seisflows.solver.Base,
"""
import os
import sys
from glob import glob

from seisflows.solver.specfem import Specfem
from seisflows.tools import unix, msg
from seisflows.tools.specfem import getpar, setpar


class Specfem2D(Specfem):
    """
    Python interface to Specfem2D.
    """
    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        super().__init__()

        self.required.par(
            "SOURCE_PREFIX", required=False, default="SOURCE", par_type=str,
            docstr="Prefix of SOURCE files in path SPECFEM_DATA. By "
                   "default, 'SOURCE' for SPECFEM2D"
        )

        self.f0 = None

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
        Returns a wildcard identifier for synthetic data based on SPECFEM2D
        file naming schema. Allows formatting dcomponent e.g.,
        when called by solver.data_filenames

        :type comp: str
        :param comp: component formatter, defaults to wildcard '?'
        :rtype: str
        :return: wildcard identifier for channels
        """
        if self.par.FORMAT.upper() == "SU":
            # return f"*.su"  # too vague but maybe for a reason? -bryant
            return f"U{comp}_file_single.su"
        elif self.par.FORMAT.upper() == "ASCII":
            return f"*.?X{comp}.sem?"

    @property
    def data_filenames(self):
        """
        Returns the filenames of all data, either by the requested components
        or by all available files in the directory.

        .. note::
            If the glob returns an  empty list, this function exits the
            workflow because filenames should  not be empty is they're being
            queried

        :rtype: list
        :return: list of data filenames
        """
        unix.cd(self.cwd)
        unix.cd(os.path.join("traces", "obs"))

        if self.par.COMPONENTS:
            filenames = []
            if self.par.FORMAT.upper() == "SU":
                for comp in self.par.COMPONENTS:
                    filenames += [self.data_wildcard.format(comp=comp.lower())]
                    # filenames += [f"U{comp.lower()}_file_single.su"]
            elif self.par.FORMAT.upper() == "ASCII":
                for comp in self.par.COMPONENTS:
                    filenames += glob(
                            self.data_wildcard.format(comp=comp.upper())
                            )
                    # filenames += glob(f"*.?X{comp.upper()}.sem?")
        else:
            filenames = glob(self.data_wildcard)

        if not filenames:
            print(msg.cli("The property solver.data_filenames, used to search "
                          "for traces in 'scratch/solver/*/traces' is empty "
                          "and should not be. Please check solver parameters: ",
                          items=[f"data_wildcard: {self.data_wildcard}"],
                          header="data filenames error", border="=")
                  )
            sys.exit(-1)

        return filenames

    @property
    def model_databases(self):
        """
        The location of model inputs and outputs as defined by SPECFEM2D
        """
        return os.path.join(self.cwd, "DATA")

    @property
    def kernel_databases(self):
        """
        The location of kernel inputs and outputs as defined by SPECFEM2D
        """
        return os.path.join(self.cwd, "OUTPUT_FILES")

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

    def setup(self):
        """
        Additional SPECFEM2D setup steps
        """
        super().setup()
        self.f0 = getpar(key="f0",
                         file=os.path.join(self.cwd, "DATA/SOURCE"))[1]

        if "MULTIPLES" in self.par:
            if self.par.MULTIPLES:
                setpar(key="absorbtop", val=".false.", file="DATA/Par_file")
            else:
                setpar(key="absorbtop", val=".true.", file="DATA/Par_file")

    def _set_model(self, model_name, model_type="gll"):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._set_model(model_name=model_name, model_type=model_type)

    def generate_data(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self.generate_data()

    def eval_func(self, path, write_residuals=True):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self.eval_func(path=path, write_residuals=write_residuals)

    def eval_grad(self, path, export_traces=True):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self.eval_grad(path=path, export_traces=export_traces)

    def _forward(self, output_path):
        """
        Calls SPECFEM2D forward solver, exports solver outputs to traces dir

        :type output_path: str
        :param output_path: path to export traces to after completion of
            simulation expected values are either 'traces/obs' for 'observation'
            data (i.e., synthetics generated by the TRUE model), or
            'traces/syn', for synthetics generated during function evaluations
        """
        unix.cd(self.cwd)

        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")

        self._call_solver(executable="bin/xmeshfem2D", output="fwd_mesher.log")
        self._call_solver(executable="bin/xspecfem2D", output="fwd_solver.log")

        if self.par.FORMAT.upper() == "SU":
            # Work around SPECFEM2D's version dependent file names
            for tag in ["d", "v", "a", "p"]:
                unix.rename(old=f"single_{tag}.su", new="single.su",
                            names=glob(os.path.join("OUTPUT_FILES", "*.su")))

        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard)),
                dst=output_path)

    def _adjoint(self):
        """
        Calls SPECFEM2D adjoint solver, creates the `SEM` folder with adjoint
        traces which is required by the adjoint solver
        """
        unix.cd(self.cwd)

        setpar(key="SIMULATION_TYPE", val="3", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".false.", file="DATA/Par_file")

        unix.rm("SEM")
        unix.ln("traces/adj", "SEM")

        # Deal with different SPECFEM2D name conventions for regular traces and
        # "adjoint" traces
        if self.par.FORMAT.upper == "SU":
            unix.rename(old=".su", new=".su.adj",
                        names=glob(os.path.join("traces", "adj", "*.su")))

        self._call_solver(executable="bin/xspecfem2D", output="adj_solver.log")

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
        """
        Specfem2D requires additional model parameters in directory to perform
        the xsmooth_sem task. This function will copy these files into the
        directory before performing the base smooth operations.

        Kwargs should match arguments of solver.base.smooth()

        .. note::
            This operation is usually run with run(single=True) so only one
            task will be performing these operations.

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
        # Redundant to 'base' class but necessary
        if not os.path.exists(input_path):
            unix.mkdir(input_path)

        unix.cd(self.cwd)
        unix.cd("DATA")

        # Copy over only the files that are required. Won't execute if no match
        files = []
        for tag in ["jacobian", "NSPEC_ibool", "x", "y", "z"]:
            files += glob(f"*_{tag}.bin")
        for src in files:
            unix.cp(src=src, dst=input_path)

        super().smooth(input_path=input_path, output_path=output_path,
                       parameters=parameters, span_h=span_h, span_v=span_v,
                       output=output)

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
        Setup utility: Creates the "adjoint traces" expected by SPECFEM.
        This is only done for the 'base' the Preprocess class.

        .. note::
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        super()._initialize_adjoint_traces()
    
        unix.cd(self.cwd)
        unix.cd(os.path.join("traces", "adj"))

        # work around SPECFEM2D's use of different name conventions for
        # regular traces and 'adjoint' traces
        if self.par.FORMAT.upper() == "SU":
            files = glob("*SU")
            unix.rename(old="_SU", new="_SU.adj", names=files)
        elif self.par.FORMAT.upper() == "ASCII":
            files = glob("*sem?")

            # Get the available extensions, which are named based on unit
            extensions = set([os.path.splitext(_)[-1] for _ in files])
            for extension in extensions:
                unix.rename(old=extension, new=".adj", names=files)

        # SPECFEM2D requires that all components exist even if ununsed
        components = ["x", "y", "z", "p"]

        if self.par.FORMAT.upper() == "SU":
            for comp in components:
                src = f"U{self.par.COMPONENTS[0]}_file_single.su.adj"
                dst = f"U{comp.lower()}s_file_single.su.adj"
                if not os.path.exists(dst):
                    unix.cp(src, dst)
        elif self.par.FORMAT.upper() == "ASCII":
            for fid in glob("*.adj"):
                net, sta, cha, ext = fid.split(".")
                for comp in components:
                    # Replace the last value in the channel with new component
                    cha_check = cha[:-1] + comp.upper()
                    fid_check = ".".join([net, sta, cha_check, ext])
                    if not os.path.exists(fid_check):
                        unix.cp(fid, fid_check)

    def _check_mesh_properties(self, path=None):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._check_mesh_properties(path=path)

    def _check_source_names(self):
        """Inherits from seisflows.solver.specfem.Specfem"""
        self._check_source_names()



