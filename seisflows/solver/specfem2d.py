#!/usr/bin/env python3
"""
This is the subclass seisflows.solver.specfem2d

This class provides utilities for the Seisflows solver interactions with
Specfem2D. It inherits all attributes from seisflows.solver.Base,
"""
import os
import sys
import logging
from glob import glob

from seisflows.tools import unix, msg
from seisflows.tools.wrappers import exists
from seisflows.config import custom_import, SeisFlowsPathsParameters
from seisflows.tools.specfem import call_solver, getpar, setpar


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
    # Class-specific logger accessed using self.logger
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
        Checks SPECFEM2D Par_file for acceptable parameters and matches with
        the internally set parameters
        """
        # Check the number of steps in the SPECFEM2D Par_file
        nt_str, nt, nt_i = getpar(key="NSTEP", file="DATA/Par_file")
        if int(nt) != PAR.NT:
            if self.taskid == 0:
                print(msg.cli(f"SPECFEM2D {nt_str}=={nt} is not equal "
                              f"SeisFlows3 PAR.NT=={PAR.NT}. Please ensure "
                              f"that these values match in both files.",
                              header="parameter match error", border="=")
                      )
                sys.exit(-1)

        dt_str, dt, dt_i = getpar(key="DT", file="DATA/Par_file")
        if float(dt) != PAR.DT:
            if self.taskid == 0:
                print(msg.cli(f"SPECFEM2D {dt_str}=={dt} is not equal "
                              f"SeisFlows3 PAR.DT=={PAR.DT}. Please ensure "
                              f"that these values match in both files.",
                              header="parameter match error", border="=")
                      )
                sys.exit(-1)

        # Check the central frequency in the SPECFEM2D SOURCE file
        f0_str, f0, f0_i = getpar(key="f0", file="DATA/SOURCE")
        if float(f0) != PAR.F0:
            if self.taskid == 0:
                print(msg.cli(f"SPECFEM2D {f0_str}=={f0} is not equal "
                              f"SeisFlows3 PAR.F0=={PAR.F0}. Please ensure "
                              f"that these values match the DATA/SOURCE file.",
                              header="parameter match error", border="=")
                      )
                sys.exit(-1)

        # Ensure that NPROC matches the MESH values
        nproc = self.mesh_properties.nproc
        if nproc != PAR.NPROC:
            if self.taskid == 0:
                print(msg.cli(f"SPECFEM2D mesh NPROC=={nproc} is not equal"
                              f"SeisFlows3 PAR.NPROC=={PAR.NPROC}. "
                              f"Please check that your mesh matches this val.",
                              header="parameter match error", border="=")
                      )
                sys.exit(-1)

        if "MULTIPLES" in PAR:
            if PAR.MULTIPLES:
                setpar(key="absorbtop", val=".false.", file="DATA/Par_file")
            else:
                setpar(key="absorbtop", val=".true.", file="DATA/Par_file")

    def generate_data(self, **model_kwargs):
        """
        Generates data using the True model, exports traces to `traces/obs`

        :param model_kwargs: keyword arguments to pass to `generate_mesh`
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.cwd)
        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")

        call_solver(PAR.MPIEXEC, "bin/xmeshfem2D", output="mesher.log")
        call_solver(PAR.MPIEXEC, "bin/xspecfem2D", output="solver.log")

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
        setpar(key="SIMULATION_TYPE", val="1", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".true.", file="DATA/Par_file")

        call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xmeshfem2D")
        call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xspecfem2D")

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
        setpar(key="SIMULATION_TYPE", val="3", file="DATA/Par_file")
        setpar(key="SAVE_FORWARD", val=".false.", file="DATA/Par_file")

        unix.rm("SEM")
        unix.ln("traces/adj", "SEM")

        # Deal with different SPECFEM2D name conventions for regular traces and
        # "adjoint" traces
        if PAR.FORMAT.upper == "SU":
            unix.rename(old=".su", new=".su.adj",
                        names=glob(os.path.join("traces", "adj", "*.su")))

        call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xmeshfem2D")
        call_solver(mpiexec=PAR.MPIEXEC, executable="bin/xspecfem2D")

    def smooth(self, input_path, **kwargs):
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
        """
        # Redundant to 'base' class but necessary
        if not exists(input_path):
            unix.mkdir(input_path)

        unix.cd(self.cwd)
        unix.cd("DATA")

        # Copy over only the files that are required. Won't execute if no match
        files = []
        for tag in ["jacobian", "NSPEC_ibool", "x", "y", "z"]:
            files += glob(f"*_{tag}.bin")
        for src in files:
            unix.cp(src=src, dst=input_path) 

        super().smooth(input_path=input_path, **kwargs)

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

        .. note:: 
            If the glob returns an  empty list, this function exits the 
            workflow because filenames should  not be empty is they're being 
            queried

        :rtype: list
        :return: list of data filenames
        """
        unix.cd(self.cwd)
        unix.cd(os.path.join("traces", "obs"))

        if PAR.COMPONENTS:
            filenames = []
            if PAR.FORMAT.upper() == "SU":
                for comp in PAR.COMPONENTS:
                    filenames += [self.data_wildcard.format(comp=comp.lower())]
                    # filenames += [f"U{comp.lower()}_file_single.su"]
            elif PAR.FORMAT.upper() == "ASCII":
                for comp in PAR.COMPONENTS:
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

    @property
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
        if PAR.FORMAT.upper() == "SU":
            # return f"*.su"  # too vague but maybe for a reason? -bryant
            return f"U{comp}_file_single.su"
        elif PAR.FORMAT.upper() == "ASCII":
            return f"*.?X{comp}.sem?"

    @property
    def source_prefix(self):
        """
        Specfem2D's preferred source prefix

        :rtype: str
        :return: source prefix
        """
        return PAR.SOURCE_PREFIX.upper()

