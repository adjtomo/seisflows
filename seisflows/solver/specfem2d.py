#!/usr/bin/env python
"""
This is the subclass seisflows.solver.specfem2d

This class provides utilities for the Seisflows solver interactions with
Specfem2D. It inherits all attributes from seisflows.solver.Base,
"""
import os
import sys
from glob import glob

from seisflows.plugins.solver.specfem2d import smooth_legacy
from seisflows.tools.seismic import getpar, setpar
from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import custom_import, SeisFlowsPathsParameters
from seisflows.tools.seismic import call_solver


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
preprocess = sys.modules['seisflows_preprocess']


class Specfem2D(custom_import("solver", "base")):
    """
    Python interface to Specfem2D. This subclass inherits functions from
    seisflows.solver.Base

    !!! See base class for method descriptions !!!
    """
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

        sf.par("F0", required=True, par_type=float,
               docstr="Dominant source frequency")

        sf.par("FORMAT", required=True, par_type=float,
               docstr="Format of synthetic waveforms used during workflow, "
                      "available options: ['ascii', 'su']")

        sf.par("SOURCE_PREFIX", required=False, default="SOURCE",
               par_type=str,
               docstr="Prefix of SOURCE files in path SPECFEM_DATA. By "
                      "default, 'SOURCE' for SPECFEM2D")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

        acceptable_formats = ["SU", "ASCII"]
        assert(PAR.FORMAT.upper() in acceptable_formats), \
            f"FORMAT must be {acceptable_formats}"

    def check_solver_parameter_files(self):
        """
        Checks solver parameters
        """
        def getpar_tryexcept(trial_list, cast, tag=""):
            """
            Re-used function to wrap getpar() in a try-except
            To allow for different SPECFEM2D Par_file version

            :type trial_list: list
            :param trial_list: list of strings to check in par file
            :type cast: str
            :param cast: cast the output results as this type
            :type tag: str
            :param tag: tag used incase error raised, for more useful message
            :rtype tuple: (str, cast)
            :return: the correct check from the trial list and corresponding val
            """
            for check in trial_list:
                try:
                    return check, getpar(check, cast=cast)
                except KeyError as e:
                    pass
            else:
                raise KeyError(f"Parameter '{tag}' not found when looking for "
                               f"{trial_list}") from e

        # Check the number of steps in the SPECFEM2D Par_file
        nt_str, nt = getpar_tryexcept(trial_list=["NSTEP", "nt"],
                                      cast=int, tag="nt")
        if nt != PAR.NT:
            if self.taskid == 0:
                print(f"WARNING: nt={nt} not equal PAR.NT={PAR.NT},"
                      f"setting PAR FILE nt={PAR.NT}")
            setpar(nt_str, PAR.NT)

        # Check the dt step discretization in the SPECFEM2D Par_file
        dt_str, dt = getpar_tryexcept(trial_list=["DT", "deltat"],
                                      cast=float, tag="dt")
        if dt != PAR.DT:
            if self.taskid == 0:
                print(f"WARNING: dt={dt} not equal PAR.DT={PAR.DT},"
                      f"setting PAR FILE dt={PAR.DT}")
            setpar(dt_str, PAR.DT)

        # Check the central frequency in the SPECFEM2D SOURCE file
        f0 = getpar("f0", file="DATA/SOURCE", cast=float)
        if f0 != PAR.F0:
            if self.taskid == 0:
                print(f"WARNING: f0={f0} not equal PAR.F0={PAR.F0},"
                      f"setting SOURCE f0={PAR.F0}")
            setpar("f0", PAR.F0, filename="DATA/SOURCE")

        # Ensure that NPROC matches the MESH values
        if self.mesh_properties.nproc != PAR.NPROC:
            if self.taskid == 0:
                print(f"Warning: "
                      f"mesh_properties.nproc={self.mesh_properties.nproc} "
                      f"not equal  PAR.NPROC={PAR.NPROC}"
                      )

        if "MULTIPLES" in PAR:
            if PAR.MULTIPLES:
                setpar("absorbtop", ".false.")
            else:
                setpar("absorbtop", ".true.")

    def generate_data(self, **model_kwargs):
        """
        Generates data using the True model, exports traces to `traces/obs`

        :param model_kwargs: keyword arguments to pass to `generate_mesh`
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.cwd)
        setpar("SIMULATION_TYPE", "1")
        setpar("SAVE_FORWARD", ".true.")

        call_solver(system.mpiexec(), "bin/xmeshfem2D", output="mesher.log")
        call_solver(system.mpiexec(), "bin/xspecfem2D", output="solver.log")

        if PAR.FORMAT.upper() == "SU":
            # Work around SPECFEM2D's version dependent file names
            for tag in ["d", "v", "a", "p"]:
                unix.rename(old=f"single_{tag}.su", new="single.su",
                            names=glob(os.path.join("OUTPUT_FILES", "*.su")))

        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard)),
                dst=os.path.join("traces", "obs"))

        if PAR.SAVETRACES:
            self.export_traces(os.path.join(PATH.OUTPUT, "traces", "obs"))

    def initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM.
        This is only done for the 'base' the Preprocess class.

        Note:
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        super().initialize_adjoint_traces()
    
        unix.cd(self.cwd)
        unix.cd(os.path.join("traces", "adj"))

        # work around SPECFEM2D's use of different name conventions for
        # regular traces and 'adjoint' traces
        if PAR.FORMAT.upper() == "SU":
            files = glob("*SU")
            unix.rename(old="_SU", new="_SU.adj", names=files)
        elif PAR.FORMAT.upper() == "ASCII":
            files = glob("*sem?")

            # Get the available extensions, which are named based on unit
            extensions = set([os.path.splitext(_)[-1] for _ in files])
            for extension in extensions:
                unix.rename(old=extension, new=".adj", names=files)

        # SPECFEM2D requires that all components exist even if ununsed
        components = ["x", "y", "z", "p"]

        if PAR.FORMAT.upper() == "SU":
            for comp in components:
                src = f"U{PAR.COMPONENTS[0]}_file_single.su.adj"
                dst = f"U{comp.lower()}s_file_single.su.adj"
                if not exists(dst):
                    unix.cp(src, dst)
        elif PAR.FORMAT.upper() == "ASCII":
            for fid in glob("*.adj"):
                net, sta, cha, ext = fid.split(".")
                for comp in components:
                    # Replace the last value in the channel with new component
                    cha_check = cha[:-1] + comp.upper()
                    fid_check = ".".join([net, sta, cha_check, ext])
                    if not exists(fid_check):
                        unix.cp(fid, fid_check)

    def generate_mesh(self, model_path, model_name, model_type='gll'):
        """
        Performs meshing with internal mesher Meshfem2D and database generation

        :type model_path: str
        :param model_path: path to the model to be used for mesh generation
        :type model_name: str
        :param model_name: name of the model to be used as identification
        :type model_type: str
        :param model_type: available model types to be passed to the Specfem3D
            Par_file. See Specfem3D Par_file for available options.
        """
        assert(exists(model_path)), f"model {model_path} does not exist"

        available_model_types = ["gll"]
        assert(model_type in available_model_types), \
            f"{model_type} not in available types {available_model_types}"

        unix.cd(self.cwd)

        # Run mesh generation
        if model_type == "gll":
            self.check_mesh_properties(model_path)

            # Copy the model files (ex: proc000023_vp.bin ...) into DATA
            src = glob(os.path.join(model_path, "*"))
            dst = self.model_databases
            unix.cp(src, dst)

        # Export the model into output folder
        if self.taskid == 0:
            self.export_model(os.path.join(PATH.OUTPUT, model_name))

    def forward(self, path='traces/syn'):
        """
        Calls SPECFEM2D forward solver, exports solver outputs to traces dir

        :type path: str
        :param path: path to export traces to after completion of simulation
        """
        setpar("SIMULATION_TYPE", "1")
        setpar("SAVE_FORWARD", ".true.")

        call_solver(mpiexec=system.mpiexec(), executable="bin/xmeshfem2D")
        call_solver(mpiexec=system.mpiexec(), executable="bin/xspecfem2D")

        if PAR.FORMAT.upper() == "SU":
            # Work around SPECFEM2D's version dependent file names
            for tag in ["d", "v", "a", "p"]:
                unix.rename(old=f"single_{tag}.su", new="single.su",
                            names=glob(os.path.join("OUTPUT_FILES", "*.su")))

        unix.mv(src=glob(os.path.join("OUTPUT_FILES", self.data_wildcard)),
                dst=path)

    def adjoint(self):
        """
        Calls SPECFEM2D adjoint solver, creates the `SEM` folder with adjoint
        traces which is required by the adjoint solver
        """
        setpar("SIMULATION_TYPE", "3")
        setpar("SAVE_FORWARD", ".false.")

        unix.rm("SEM")
        unix.ln("traces/adj", "SEM")

        # Deal with different SPECFEM2D name conventions for regular traces and
        # "adjoint" traces
        if PAR.FORMAT.upper == "SU":
            unix.rename(old=".su", new=".su.adj",
                        names=glob(os.path.join("traces", "adj", "*.su")))

        call_solver(mpiexec=system.mpiexec(), executable="bin/xmeshfem2D")
        call_solver(mpiexec=system.mpiexec(), executable="bin/xspecfem2D")

    def import_model(self, path):
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

    def export_model(self, path):
        """
        File transfer utility to move a SPEFEM2D model from the DATA directory
        to an external path location

        :type path: str
        :param path: path to export the SPECFEM2D model
        :return:
        """
        unix.mkdir(path)
        unix.cp(src=glob(os.path.join(self.cwd, "DATA", "*.bin")),
                dst=path)

    @property
    def data_filenames(self):
        """
        Returns the filenames of all data, either by the requested components
        or by all available files in the directory.

        :rtype: list
        :return: list of data filenames
        """
        unix.cd(self.cwd)
        unix.cd(os.path.join("traces", "obs"))

        if PAR.COMPONENTS:
            filenames = []
            if PAR.FORMAT.upper() == "SU":
                for comp in PAR.COMPONENTS:
                    filenames += [f"U{comp.lower()}_file_single.su"]
            elif PAR.FORMAT.upper() == "ASCII":
                for comp in PAR.COMPONENTS:
                    filenames += glob(f"*.?X{comp.upper()}.sem?")
            return filenames
        else:
            return glob(self.data_wildcard)

    @property
    def model_databases(self):
        """
        The location of databases for kernel outputs
        """
        return os.path.join(self.cwd, "DATA")

    @property
    def kernel_databases(self):
        """
        The location of databases for model outputs
        """
        return os.path.join(self.cwd, "OUTPUT_FILES")

    @property
    def data_wildcard(self):
        """
        Returns a wildcard identifier for synthetic data

        :rtype: str
        :return: wildcard identifier for channels
        """
        if PAR.FORMAT.upper() == "SU":
            # return f"*.su"  # too vague but maybe for a reason? -bryant
            return f"U?_file_single.su"
        elif PAR.FORMAT.upper() == "ASCII":
            return f"*.?X?.sem?"

    @property
    def source_prefix(self):
        """
        Specfem2D's preferred source prefix

        :rtype: str
        :return: source prefix
        """
        return PAR.SOURCE_PREFIX.upper()

