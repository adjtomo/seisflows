#!/usr/bin/env python3
"""
A command line tool for using and manipulating SeisFlows.
The main entry point to the SeisFlows package, this command line tool
facilitates interface with the underlying SeisFlows package.

.. rubric::
    $ seisflows -h  # runs the help command to investigate package features

.. note::
    To add new functions to the seisflows command line tool, you must:
    - Write a new function within the SeisFlows class
    - Add a new subparser with optional arguments to sfparser()
    - Add subparser to subparser dict at the end of sfparser()
    In-function import statements are used to reduce call-time for simpler fx's
"""
import os
import sys
import argparse
from glob import glob
from inspect import getmro
from seisflows import logger, ROOT_DIR, NAMES
from seisflows.tools import unix, msg
from seisflows.tools.config import load_yaml, custom_import, import_seisflows
from seisflows.tools.specfem import (getpar, setpar, getpar_vel_model,
                                     setpar_vel_model)



def sfparser():
    """
    An command-line argument parser which allows for intuitive exploration of
    the available functions.

    Gets User defined arguments or assign defaults. Makes use of subparsers to
    get individual help statements for each of the main functions.

    .. rubric::
        $ seisflows {main arg} {optional sub arg}

    :rtype: argparse.ArgumentParser()
    :return: User defined or default arguments
    """
    class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
        """
        Override the help statement to NOT print out available subcommands for a
        cleaner UI when calling this CLI tool.

        https://stackoverflow.com/questions/13423540/
                              argparse-subparser-hide-metavar-in-command-listing
        """
        def _format_action(self, action):
            parts = super()._format_action(action)
            if action.nargs == argparse.PARSER:
                parts = "\n".join(parts.split("\n")[1:])
            return parts

    # Initiate the argument parser with a nicely formatted ASCII descriptor
    parser = argparse.ArgumentParser(
        formatter_class=SubcommandHelpFormatter,
        description=f"{'='*80}\n\n"
                    f"{'SeisFlows: Waveform Inversion Package':^80}\n\n"
                    f"{'='*80}",
        epilog="'seisflows [command] -h' for more detailed descriptions "
               "of each command.",
    )

    # Optional parameters
    parser.add_argument("-w", "--workdir", nargs="?", default=os.getcwd(),
                        help="The SeisFlows working directory, default: cwd")
    parser.add_argument("-p", "--parameter_file", nargs="?",
                        default="parameters.yaml",
                        help=f"Parameters file, default: 'parameters.yaml'")

    # Initiate a sub parser to provide nested help functions and sub commands
    subparser = parser.add_subparsers(
        title="command",
        description="Available SeisFlows arguments and their intended usages",
        dest="command",
    )
    # The following subparsers constitute the available SeisFlows commands
    # and each refers to a function within the SeisFlows class.
    # =========================================================================
    setup = subparser.add_parser(
        "setup", help="Setup working directory from scratch",
        description="""In the specified working directory, copy template 
        parameter file containing only module choices, and symlink source code 
        for both the base and super repositories for easy edit access. If a 
        parameter file matching the provided name exists in the working 
        directory, a prompt will appear asking the user if they want to 
        overwrite."""
    )
    setup.add_argument("-f", "--force", action="store_true",
                       help="automatically overwrites existing parameter file")
    # =========================================================================
    configure = subparser.add_parser(
        "configure", help="Fill parameter file with defaults",
        description="""SeisFlows parameter files will vary depending on 
        chosen modules and their respective required parameters. This function 
        will dynamically traverse the source code and generate a template 
        parameter file based on module choices. The resulting file incldues 
        docstrings and type hints for each parameter. Optional parameters will 
        be set with default values and required parameters and paths will be 
        marked appropriately. Required parameters must be set before a workflow
        can be submitted."""
    )
    configure.add_argument("-a", "--absolute_paths", action="store_true",
                           help="Set default paths relative to cwd")
    # =========================================================================
    swap = subparser.add_parser(
        "swap", help="Swap module parameters in an existing parameter file",
        description="""During workflow development, it may be necessary to swap
        between different sub-modules (e.g., system.workstation -> 
        system.cluster). However this would typically involving re-generating
        and re-filling a parameter file. The 'swap' function makes it easier
        to swap parameters between modules.
        """
    )
    swap.add_argument("module", nargs="?", help="Module name to swap")
    swap.add_argument("classname", nargs="?", help="Classname to swap to")
    # =========================================================================
    submit = subparser.add_parser(
        "submit", help="Submit initial workflow to system",
        description="""The main SeisFlows execution command. Submit a SeisFlows 
        workflow to the chosen system, equal to executing 
        seisflows.workflow.main(). This function will create and fill the 
        working directory with required paths, perform path and parameter 
        error checking, and establish the active working environment before
        executing the workflow."""
    )
    submit.add_argument("-s", "--stop_after", default=None, type=str,
                        help="Optional override of the 'STOP_AFTER' parameter")
    # =========================================================================
    resume = subparser.add_parser(
        "resume", help="Re-submit previous workflow to system",
        description="""Resume a previously submitted workflow. Used when 
        an active environment exists in the working directory, and must be 
        submitted to the system again."""
    )
    resume.add_argument("-r", "--resume_from", default=None, type=str,
                        help="Optional override of the 'RESUME_FROM' parameter")
    resume.add_argument("-s", "--stop_after", default=None, type=str,
                        help="Optional override of the 'STOP_AFTER' parameter")
    # =========================================================================
    restart = subparser.add_parser(
        "restart", help="Remove current environment and submit new workflow",
        description="""Akin to running seisflows clean; seisflows submit. 
        Restarts the workflow by removing the current state and submitting a 
        fresh workflow."""
    )
    restart.add_argument("-f", "--force", action="store_true",
                         help="Skip the clean warning check statement")
    # =========================================================================
    clean = subparser.add_parser(
        "clean", help="Remove files relating to an active working environment",
        description="""Delete all SeisFlows related files in the working 
        directory, except for the parameter file."""
    )
    clean.add_argument("-f", "--force", action="store_true",
                       help="Skip the warning check that precedes the clean "
                       "function")
    # =========================================================================
    par = subparser.add_parser(
        "par", help="View and edit SeisFlows parameter file",
        description="""Directly edit values in the parameter file by providing
        the parameter and corresponding value. If no value is provided, will 
        simply print out the current value of the given parameter. Works also
        with path names."""
    )
    par.add_argument("parameter", nargs="?", help="Parameter to edit or view, "
                     "(case independent).")
    par.add_argument("value", nargs="?", default=None,
                     help="Optional value to set parameter to. If not given, "
                     "will print out current parameter. If given, will replace "
                     "current parameter with new value. Set as 'null' "
                     "for NoneType and set '' for empty string")
    par.add_argument("-p", "--skip_print", action="store_true", default=False,
                     help="Skip the print statement which is typically "
                          "sent to stdout after changing parameters.")
    # =========================================================================
    sempar = subparser.add_parser(
        "sempar", help="View and edit SPECFEM parameter file",
        description="""Directly edit values in the SPECFEM parameter file by 
        providing the parameter and corresponding value. If no value is 
        provided, will simply print out the current value of the given 
        parameter. Works also with path names."""
    )
    sempar.add_argument("parameter", nargs="?", help="Parameter to edit or "
                        "view (case independent)")
    sempar.add_argument("value", nargs="?", default=None,
                        help="Optional value to set parameter to.")
    sempar.add_argument("-P", "--par_file", nargs="?", default="Par_file",
                        help="Parameter file")
    sempar.add_argument("-p", "--skip_print", action="store_true",
                        default=False,
                        help="Skip the print statement which is typically "
                             "sent to stdout after changing parameters.")

    # =========================================================================
    check = subparser.add_parser(
        "check",  formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Run check functions to ensure that the provided parameter "
                    "file has been set correctly"
    )
    # =========================================================================
    init = subparser.add_parser(
        "init", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Run check and setup functions to generate a SeisFlows "
                    "working directory")
    # =========================================================================
    plot2d = subparser.add_parser(
        "plot2d", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Plots model/kernels/gradient files located in the output/
        directory. ONLY available for SPECFEM2D models.""",
        help="Plot 2D figures of models/kernels/gradients.")

    plot2d.add_argument("name", type=str, nargs="?",
                        help="Name of directory in the output/ directory")
    plot2d.add_argument("parameter", type=str, nargs="?",
                        help="Name of parameter to plot from `name`. E.g., 'vs', "
                             "'vp' etc.")
    plot2d.add_argument("-c", "--cmap", type=str, nargs="?",
                        help="colormap to be passed to PyPlot")
    plot2d.add_argument("-s", "--savefig", type=str, nargs="?", default=None,
                        help="optional name and path to save figure")
    # =========================================================================
    plotst = subparser.add_parser(
        "plotst", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Plots waveforms output by the solver. Uses ObsPy's 
Stream.plot() function under the hood. Example call would be 
`seisflows plotst scratch/solver/mainsolver/traces/syn/*`
        """)

    plotst.add_argument("fids", type=str, nargs="*",
                        help="File IDs to be passed to plotting. Wildcards "
                             "acceptable")
    plotst.add_argument("--data_format", type=str, nargs="?", default="ASCII",
                        help="Data format of the files. Must match file type "
                             "that SeisFlows can read. These include:"
                             "['SU', 'ASCII']. Defaults to 'ASCII'. See "
                             "SeisFlows.preprocess.default.read() for "
                             "all options.")
    plotst.add_argument("-s", "--savefig", type=str, nargs="?", default=None,
                        help="optional name and path to save figure")
    # =========================================================================
    print_ = subparser.add_parser(
        "print", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Print information related to an active environment

    modules       List available module names for all available packages
    flow          Print out the workflow.main() flow arguments
    inherit       Track inheritance chain for all modules, determine method 
                  ownership for a given function. 
                  seisflows print inherit {optional module} {optional function}
                  e.g., seisflows inherit workflow main
                    """,
        help="Print information related to an active environment")

    print_.add_argument("choice", type=str, nargs="?",
                        help="Parameter, state, or value to check")
    print_.add_argument("args", type=str, nargs="*",
                        help="Generic arguments passed to check functions")
    # =========================================================================
    reset = subparser.add_parser(
        "reset", formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Reset modules within an active state", description="""
Occasionally the machinery of a given module must be reset within an active 
working state before the workflow can be resumed

    line_search     Reset line search, step count returns to 1
                         """)
    reset.add_argument("choice", type=str, nargs="?", default=None,
                       help="Choice of module/component to reset")
    reset.add_argument("args", type=str, nargs="*",
                       help="Generic arguments passed to reset functions")
    # =========================================================================
    subparser.add_parser(
        "debug", help="Start interactive debug environment",
        description="""Starts an IPython debugging environment and loads an
        active SeisFlows working state, as well as distributing the SeisFlows
        module namespace. Allows exploration of the active state, as well as
        manually control of the workflow. Useful for recovery from unexpected
        workflow crashes. State changes will not be saved automatically. Type
        'workflow.checkpoint()' in the debug environment to save any changes
        made during debugging.
        """)
    # =========================================================================
    examples = subparser.add_parser(
        "examples", help="Look at and run pre-configured example problems",
        description="""Lists out available example problems and allows the
        user to run example problems directly from the command line. Some 
        example problems may have pre-run prompts mainly involving the
        numerical solver
        """
    )
    examples.add_argument("method", type=str, nargs="?", default=None,
                          help="Method for running the example problem. If not"
                               "provided, simply prints out the list of "
                               "available example problems. If given as an "
                               "integer value, will print out the help message "
                               "for the given example. If 'run', will run the "
                               "example. If 'setup' will simply setup the "
                               "example working directory but will not execute "
                               "`seisflows submit`")
    examples.add_argument("choice", type=int,  nargs="?", default=None,
                          help="If `method` in ['setup', 'run'], integer"
                               "value corresponding to the given example "
                               "problem which can listed using `seisflows "
                               "examples`")
    examples.add_argument("-r", "--specfem2d_repo", type=str,  nargs="?",
                          default=None,
                          help="path to the SPECFEM2D directory which should "
                               "contain binary executables. If not given, "
                               "assumes directory is called 'specfem2d/' in "
                               "the current working directory. If that dir "
                               "is not found, SPECFEM2D will be downloaded, "
                               "configured and compiled automatically in the "
                               "current working directory.")
    examples.add_argument("--nsta", type=int, nargs="?", default=None,
                          help="User-defined number of stations to use for "
                               "the example problem (1 <= NSTA <= 131). If "
                               "not given, each example has its own default.")
    examples.add_argument("--ntask", type=int, nargs="?", default=None,
                          help="User-defined number of events to use for "
                               "the example problem (1 <= NTASK <= 25). If "
                               "not given, each example has its own default.")
    examples.add_argument("--niter", type=int, nargs="?", default=None,
                          help="User-defined number of iterations to run for "
                               "the example problem (1 <= NITER <= inf). If "
                               "not given, each example has its own default.")
    examples.add_argument("--event_id", type=int, nargs="?", default=None,
                          help="Allow User to choose a specific event ID from "
                               "the Tape 2007 example (1 <= EVENT_ID <= 25). "
                               "If not used, example will default to choosing "
                               "sequential from 1 to NTASK")
    examples.add_argument("--with_mpi", action="store_true", default=False,
                          help="Run Solver with MPI using MPI exectuable "
                               "`mpiexec` (defaults to 'mpirun'). Defaults to "
                               "False, in which case executables are run "
                               "in serial")
    examples.add_argument("--mpiexec", type=str, nargs="?", default="mpirun",
                          help="Only for Example(s) 4. MPI executable to use "
                               "when running SPECFEM2D with MPI. Defaults to "
                               "'mpirun'")
    examples.add_argument("--nproc", type=int, nargs="?", default=None,
                          help="Number of processors to use for MPI runs. "
                               "Default values chosen for each example problem."
                          )
    # =========================================================================
    # Defines all arguments/functions that expect a sub-argument
    subparser_dict = {"check": check, "par": par,
                      "sempar": sempar, "clean": clean, "plot2d": plot2d,
                      "restart": restart, "print": print_, "reset": reset,
                      "examples": examples, "swap": swap}
    if parser.parse_args().command in subparser_dict:
        return parser, subparser_dict[parser.parse_args().command]
    else:
        return parser, None


class SeisFlows:
    """
    The main entry point to the SeisFlows package, to be interacted with
    through the command line. This class is responsible for:
        1) setting up or re-creating a SeisFlows working enviornment,
        2) (re-)submitting workflows to the system,
        3) inspecting, manipulating or viewing a live working environment via
            command line arguments.

    .. rubric::
        $ seisflows -h

    .. note::
        Almost every modules requires loading of other modules, i.e. to run
        any checks we must load the entire SeisFlows environment, which is slow
        but provides the most flexibility when accessing internal information
    """
    def __init__(self, workdir=None, parameter_file=None):
        """
        Parse user-defined arguments and establish internal parameters used to
        control which functions execute and how. Instance must be called to
        execute internal functions

        .. note::
            Normally the `workdir` and `parameter_file` paramters should be
            set through the argparser (i.e., command line) but we allow
            it to be set in a Python environment through __init__ incase there
            is no access to the argparser through the command line (e.g., when
            testing functionality using pytest)

        :type workdir: str
        :param workdir: the working directory to initiate a SeisFlows workflow.
        :type parameter_file: str
        :param parameter_file: full path to the parameter file
        """
        self._parser, self._subparser = sfparser()
        self._paths = None
        self._parameters = None
        self._args = self._parser.parse_args()
        if workdir is not None:
            self._args.workdir = workdir
        if parameter_file is not None:
            self._args.parameter_file = parameter_file

    def __call__(self, command=None, **kwargs):
        """
        When called, SeisFlows will execute one of its internal functions

        .. rubric::
            # From the command line
            $ seisflows {command} {optional subcommand}

            # From inside a Python environment
            > from seisflows.scripts.seisflows import SeisFlows
            > sf = SeisFlows()
            > sf("{command}", {optional subcommand}={value})

            # Example
            $ seisflows par linesearch

        :type command: str
        :param command: If not None, allows controlling this class from inside
            a Python environment. If sub-commands are required, these are
            inserted using the kwargs.
            Usually not required unless writing tests or scripting SF in Python
        :type return_self: bool
        :param return_self: if True, do not execute a command, which init
            usually does, but return the SeisFlows class itself. This is used
            just for testing purposes
        :return:
        """
        if command is not None:
            # This allows running SeisFlows() from inside a Python environment
            # mostly used for testing purposes but can also be used for scripts
            kwargs = {**kwargs, **vars(self._args)}  # include argparse defaults
            getattr(self, command)(**kwargs)
        else:
            # This is the main command-line functionality of the class
            # Print out the help statement if no command is given
            if len(sys.argv) == 1:
                self._parser.print_help()
                sys.exit(0)

            # Call the given function based on the user-defined name.
            # Throw in all arguments as kwargs and let the function sort it out
            getattr(self, self._args.command)(**vars(self._args))

    @property
    def _public_methods(self):
        """
        Return a list of all public methods within this class.

        .. warning::
            Only methods that can be called via the command line should be
            public, all other methods and attributes should be private.
        """
        return [_ for _ in dir(self) if not _.startswith("_")]

    def setup(self, force=False, **kwargs):
        """
        Initiate a SeisFlows working directory from scratch by establishing a
        template parameter file.

        .. note::
            Future working directory setup functions can be placed here

        :type force: bool
        :param force: flag to force parameter file overwriting
        """
        par_file = os.path.join(self._args.workdir, self._args.parameter_file)
        if os.path.exists(par_file):
            if force:
                check = "y"
            else:
                check = input(
                    msg.cli(f"Existing parameter file "
                            f"({self._args.parameter_file}) found. Do you "
                            f"wish to overwrite with a blank file? (y/[n])"
                            ))
            if check != "y":
                sys.exit(0)

        with open(par_file, "w") as f:
            f.write(msg.base_parameter_file)
        print(msg.cli(f"created parameter file: {self._args.parameter_file}"))

    def configure(self, absolute_paths=False, **kwargs):
        """
        Dynamically generate the parameter file by writing out docstrings and
        default values for each of the SeisFlows module parameters.
        This function writes files manually, consistent with the .yaml format.

        .. note::
            This function relies on docstrings being formatted the same way
            throughout the package. Note the trailing '***' character at the end
            of the docstring. This is required for `configure` to know where one
            docstring ends and another beings. The formatting looks like:

            Title
            -----
            some description

            Parameters
            ----------
            :type a: int
            :param a: parameter a

            Paths
            -----
            :type path_a: str
            :param path_a: path for a
            ***

        :type absolute_paths: bool
        :param absolute_paths: if True, expand pathnames to absolute paths,
            else if False, use path names relative to the working directory.
            Defaults to False, uses relative paths.
        """
        from traceback import format_exc

        print("configuring SeisFlows parameter file")

        def split_module_docstring(mod, idx):
            """
            Since our docstrings are concatenated, we need to break them
            and remove the path docstrings, those come later.

            :type idx: int
            :param idx: 0 returns parameter docstrings, 1 returns path docstring
            """
            docstring = mod.__doc__.replace("\n", "\n#")
            docssplit = docstring.split("***\n#")
            docfinal = "".join([_.split("Paths\n#    -----\n#")[idx] for _ in
                                docssplit])
            return docfinal

        # Load in a barebones parameter file and instantiate specific classes
        parameters = load_yaml(os.path.join(self._args.workdir,
                                            self._args.parameter_file))
        modules = [custom_import(name, parameters[name])() for name in NAMES]

        # If writing to parameter file fails for any reason, the file will be
        # mangled, create a temporary copy that can be re-instated upon failure
        temp_par_file = f".{self._args.parameter_file}"
        unix.cp(self._args.parameter_file, temp_par_file)

        try:
            written, path_docstrings = [], []
            f = open(self._args.parameter_file, "a")
            # Write all module parameters and corresponding docstrings
            for module in modules:
                if not module:
                    continue
                docstring = split_module_docstring(module, 0)
                f.write(f"# {'=' * 77}\n#{docstring}\n# {'=' * 77}\n")
                # Write the parameters, make sure to not have the same one twice
                for key, val in vars(module).items():
                    # Skip already written, hidden vars, and paths
                    if (key in written) or key.startswith("_") or key == "path":
                        continue
                    # YAML wants NoneType to be 'null'
                    if val is None:
                        val = "null"
                    f.write(f"{key}: {val}\n")
                    written.append(key)

            # Write docstrings for publically accesible path structure
            f.write(f"# {'=' * 77}\n")
            f.write("#\n")
            f.write("#\t Paths\n")
            f.write("#\t -----\n")
            for module in modules:
                if not module:
                    continue
                docstring = split_module_docstring(module, -1)
                f.write(f"#{docstring}\n")
            f.write(f"# {'=' * 77}\n")

            # Write values for publically accessible path structure
            written = []
            for module in modules:
                if not module:
                    continue
                for key, val in module.path.items():
                    # '_key' means hidden path so don't include in par file
                    if key in written or key.startswith("_"):
                        continue
                    if val is None:
                        val = "null"
                    if absolute_paths:
                        val = os.path.abspath(val)
                    f.write(f"path_{key}: {val}\n")
                    written.append(key)
        except Exception:
            unix.rm(self._args.parameter_file)
            unix.cp(temp_par_file, self._args.parameter_file)
            logger.critical(
                msg.cli(text="seisflows configure traceback", header="error")
            )
            print(format_exc())
            raise
        else:
            unix.rm(temp_par_file)

        f.close()

    def swap(self, module, classname, **kwargs):
        """
        Swap the parameters of an existing parameter file with a new module.
        Useful for changing out parameters without having to re-make a
        parameter file from scratch. e.g., to swap systems from a workstation
        to a cluster

        TODO figure out how to match paths too

        .. rubric::
            $ seisflows swap system slurm
        """
        if module not in NAMES:
            print(msg.cli(text=f"{module} does not match {NAMES}",
                          header="error"))
            sys.exit(-1)

        # Load in old parameter file and then move it to a hidden file
        ogpars = load_yaml(self._args.parameter_file)
        unix.mv(self._args.parameter_file, f"_{self._args.parameter_file}")
        try:
            # Create a new parameter file with updated module
            self.setup(force=True)
            for name in NAMES:
                setpar(key=name, val=ogpars[name],
                       file=self._args.parameter_file, delim=":")

            # Overwrite with new parameters
            setpar(key=module, val=classname, file=self._args.parameter_file,
                   delim=":")
            self.configure()

            ogpars.pop(module)  # don't edit the parameter were changing
            for key, val in ogpars.items():
                if val is None:
                    val = "null"
                try:
                    setpar(key=key, val=val, file=self._args.parameter_file,
                           delim=":")
                except KeyError:
                    continue
        # Replace parameter file if any errors happen
        except Exception as e:
            unix.mv(f"_{self._args.parameter_file}", self._args.parameter_file)
            raise(e)
        finally:
            unix.rm(f"_{self._args.parameter_file}")

    def check(self, **kwargs):
        """
        Run check() functions for a given parameter file and each of the
        SeisFlows modules, ensuring that parameters are acceptable for the
        given set of user-defined parameters
        """
        unix.mkdir(self._args.workdir)
        unix.cd(self._args.workdir)

        workflow = import_seisflows(workdir=self._args.workdir,
                                    parameter_file=self._args.parameter_file)
        try:
            workflow.check()
        except AssertionError as e:
            print(msg.cli(str(e), border="=", header="parameter errror"))

    def init(self, **kwargs):
        """
        Run check() + setup() functions for a given parameter file and each of
        the SeisFlows modules, ensuring that parameters are acceptable for the
        given set of user-defined parameters and running setup procedure
        which may create directories and perform some file management.
        """
        unix.mkdir(self._args.workdir)
        unix.cd(self._args.workdir)

        workflow = import_seisflows(workdir=self._args.workdir,
                                    parameter_file=self._args.parameter_file)
        try:
            workflow.check()
            workflow.setup()
        except AssertionError as e:
            print(msg.cli(str(e), border="=", header="parameter errror"))

    def submit(self, **kwargs):
        """
        Main SeisFlows execution command. Submit the SeisFlows workflow to
        the chosen system, and execute seisflows.workflow.main(). Will create
        the working directory and any required paths and ensure that all
        required paths exist.
        """
        unix.mkdir(self._args.workdir)
        unix.cd(self._args.workdir)

        parameters = load_yaml(self._args.parameter_file)
        system = custom_import("system", parameters.system)(**parameters)
        system.submit(workdir=self._args.workdir,
                      parameter_file=self._args.parameter_file)

    def clean(self, force=False, **kwargs):
        """
        Clean the SeisFlows working directory except for the parameter file.

        :type force: bool
        :param force: ignore the warning check that precedes the clean()
            function, useful if you don't want any input messages popping up
        """

        # Check if the filepaths exist
        if not os.path.exists(self._args.parameter_file):
            print(msg.cli(f"SeisFlows parameter file not found: "
                          f"'{self._args.parameter_file}'. Run 'seisflows "
                          f"setup' to create a new parameter file.")
                  )
            sys.exit(-1)

        if force:
            check = "y"
        else:
            check = input(msg.cli("This will remove all workflow objects "
                                  "leaving only the parameter file. Are "
                                  "you sure you want to continue? "
                                  "(y/[n])", header="clean", border="="))

        if check == "y":
            pars = load_yaml(self._args.parameter_file)
            for name in ["scratch", "output", "log_files", "state_file", 
                         "output_log"]:
                path = f"path_{name}"
                if path in pars:
                    unix.rm(pars[path])

    def restart(self, force=False, **kwargs):
        """
        Restart simply means clean the workding dir and submit a new workflow.
        """
        self.clean(force=force)
        self.submit()

    def debug(self, **kwargs):
        """
        Initiate an IPython debugging environment to explore the currently
        active SeisFlows environment. Reloads the system modules in an
        interactive environment allowing exploration of the package space.
        Does not allow stepping through of code (not a breakpoint).
        """
        from IPython import embed

        workflow = import_seisflows(workdir=self._args.workdir,
                                    parameter_file=self._args.parameter_file)

        # Break out sub-modules and parameters so they're more easily accesible
        parameters = load_yaml(self._args.parameter_file)
        system, solver, preprocess, optimize = workflow._modules.values()

        print("Loaded SeisFlows Modules:")
        for module in [workflow, system, solver, preprocess, optimize]:
            print(f"{module.__class__}")

        print(msg.cli("SeisFlows's debug mode is an embedded IPython "
                      "environment. All modules are loaded by default. "
                      "Any changes made here will not be saved unless "
                      "you explicitely run: 'workflow.checkpoint()'",
                      header="debug", border="="))

        embed(colors="Neutral")

    def sempar(self, parameter, value=None, skip_print=False,
               par_file="Par_file", **kwargs):
        """
        check or set parameters in the SPECFEM parameter file.
        By default assumes the SPECFEM parameter file is called 'Par_file'
        But this can be overwritten by using the '-p' flag.

        usage

            seisflows sempar [parameter] [value]

            to check the parameter 'nproc' from the command line:

                seisflows sempar nstep

            to set the parameter 'model' to 'GLL':

                seisflows sempar model GLL

            to check the values of a velocity model (SPECFEM2D)

                seisflows sempar velocity_model

            to edit the values of a velocity model (SPECFEM2D)

                seisflows sempar velocity_model \
                    "1 1 2600.d0 5800.d0 3500.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0"

                OR for a two-layered model

                seisflows sempar velocity_model \
                "1 1 2600.d0 5800.d0 3500.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0 + \
                 1 1 2600.d0 5800.d0 3500.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0"

                .. note::
                    For multi-layered models, the delimiter " + " is important,
                    you must have the whitespace on either side else the
                    function won't recognize these are separate layers.

        :type parameter: str
        :param parameter: parameter to check in parameter file. case insensitive
        :type value: str
        :param value: value to set for parameter. if none, will simply print out
            the current parameter value. to set as nonetype, set to 'null'
            SPECFEM2D: if set to 'velocity_model' allows the user to set and
            edit the velocity model defined in the SPECMFE2D Par_file. Not a
            very smart capability, likely easier to do this manually.
        :type par_file: str
        :param par_file: name of the SPECFEM parameter file, defaults: Par_file
        :type skip_print: bool
        :param skip_print: skip the print statement which is typically sent
            to stdout after changing parameters.
        """
        if not os.path.exists(par_file):
            sys.exit(f"\n\tparameter file '{par_file}' does not exist\n")
        if parameter is None:
            self._subparser.print_help()
            sys.exit(0)

        # SPECFEM parameter file has both upper and lower case parameters,
        # force upper just for string checking
        parameter = parameter.upper()
        items = []  # for stdout printing

        # Use the specfem tool to grab related information
        # Special case where the velocity model in SPECFEM2D doesn't isnt
        # formatted the same as the rest of the file
        if parameter == "VELOCITY_MODEL":
            key = parameter
            items = getpar_vel_model(file=par_file)
            cur_val = ""
        else:
            try:
                key, cur_val, _ = getpar(key=parameter, file=par_file,
                                         delim="=")
            except KeyError:
                print(msg.cli(f"'{parameter}' not found in {par_file}"))
                return

        # Option 1: Simply print out the value of the given parameter
        if value is None:
            if not skip_print:
                print(msg.cli(f"{key}: {cur_val}", items=items))
        # Option 2: Replace value with user-defined input
        else:
            if parameter == "VELOCITY_MODEL":
                setpar_vel_model(file=par_file, model=value.split("+"))
                if not skip_print:
                    items.append("->")
                    items += getpar_vel_model(file=par_file)
                    print(msg.cli(f"{key}:", items=items))
            else:
                setpar(key=parameter, val=value, file=par_file, delim="=")
                if not skip_print:
                    print(msg.cli(f"{key}: {cur_val} -> {value}"))

    def par(self, parameter, value=None, skip_print=False, **kwargs):
        """
        Check or set parameters in the seisflows parameter file.

        USAGE

            seisflows par [parameter] [value]

            to check the parameter 'NPROC' from the command line:

                seisflows par nproc

            to set the parameter 'BEGIN' to 2:

                seisflows par begin 2

            to change the scratch path to the current working directory, don't
            print to stdout:

                seisflows par scratch ./ -p

        :type parameter: str
        :param parameter: parameter to check in parameter file. case insensitive
        :type value: str
        :param value: value to set for parameter. if None, will simply print out
            the current parameter value. to set as nonetype, set to 'null'
        :type skip_print: bool
        :param skip_print: skip the print statement which is typically sent
            to stdout after changing parameters.
        """
        from warnings import warn

        if not os.path.exists(self._args.parameter_file):
            sys.exit(f"\n\tparameter file '{self._args.parameter_file}' "
                     f"does not exist\n")

        if parameter is None:
            self._subparser.print_help()
            sys.exit(0)

        # SeisFlows parameter file dictates upper-case parameters
        parameter = parameter.upper()
        if isinstance(value, str) and value.lower() == "none":
            warn("to set values NoneType, use 'null' not 'none'", UserWarning)

        # Use the specfem tool to grab related information
        try:
            key, cur_val, i = getpar(key=parameter,
                                     file=self._args.parameter_file,
                                     delim=":")
        except KeyError:
            print(msg.cli(f"'{parameter}' not found in "
                          f"{self._args.parameter_file}"))
            return

        # Option 1: Simply print out the value of the given parameter
        if value is None:
            if not skip_print:
                print(msg.cli(f"{key}: {cur_val}"))
        # Option 2: Replace value with user-defined input
        else:
            setpar(key=parameter, val=value, file=self._args.parameter_file,
                   delim=":")
            if not skip_print:
                print(msg.cli(f"{key}: {cur_val} -> {value}"))

    def examples(self, method=None, choice=None, **kwargs):
        """
        List or run a SeisFlows example problems

        USAGE

            seisflows examples [run] [choice]

            To list available examples:

                seisflows examples

            To run a specific example (this is the same as 'python example.py')

                seisflows examples run 1

        :type method: bool
        :param method: if True, run an example of choice `choice`
        :type choice: str
        :param choice: The choice of example, must match the given tag or file
            name that is assigned to it
        """
        # e.g., $ seisflows examples
        if method is None:
            self._print_examples()
            sys.exit(0)
        # e.g., $ seisflows examples 1
        elif method and (choice is None):
            try:
                choice = int(method)
            except ValueError:
                print(f"`method` argument must be 'run', 'setup' or an integer "
                      f"value corresponding to one of the available examples")
                sys.exit(0)

        # Allow the examples to be dynamically recovered based on user choice
        if choice == 1:
            from seisflows.examples.ex1_homogeneous_halfspace \
                import SFExample2D as Example
        elif choice == 2:
            from seisflows.examples.ex2_hh_w_pyatoa \
                import SFPyatoaEx2D as Example
        elif choice == 3:
            from seisflows.examples.ex3_fwd_solver import SFFwdEx2D as Example
        else:
            print(f"no SeisFlows example matching given number: {choice}")
            sys.exit(0)

        # Run or setup example, or just print system dialogue
        example = Example(method=method, **kwargs)
        example.print_dialogue()

        # e.g., $ seisflows examples run 1
        if method in ["setup", "run"]:
            example.main()

    @staticmethod
    def _print_examples():
        """Simply print a list of available examples which match the format
        ex_?*.py"""
        # Gather all the available examples in the repository
        examples_dir = os.path.join(ROOT_DIR, "examples")
        examples_list = []
        example_names = sorted(glob(os.path.join(examples_dir, "ex*.py")))

        for i, fid in enumerate(example_names):
            example_name = os.path.splitext(os.path.basename(fid))[0]
            examples_list.append((i+1, example_name, fid))

        items = [f"{j}: {exname}" for j, exname, fid in examples_list]
        print(msg.cli(
            "Example options where <name_or_idx> is either the example name "
            "or corresponding index, provided below.",
            items=[
                "'seisflows examples <name_or_idx>': print example description",
                "'seisflows examples setup <name_or_idx>': setup example but "
                "don't run workflow 'seisflows examples run <name_or_idx>': "
                "setup and run example"
            ],
            header="seisflows examples"
        ))
        print(msg.cli(items=items))

    def print(self, choice=None, **kwargs):
        """
        Print information relating to an active SeisFlows environment.
        Type 'seisflows check --help' for a detailed help message.

        :type choice: str
        :param choice: underlying sub-function to choose
        """
        acceptable_args = {"modules": self._print_modules,
                           "tasks": self._print_tasks,
                           "inherit": self._print_inheritance}

        # Ensure that help message is thrown for empty commands
        if choice not in acceptable_args.keys():
            self._subparser.print_help()
            sys.exit(0)

        acceptable_args[choice](*self._args.args, **kwargs)

    def plotst(self, fids, data_format="ASCII", savefig=None, **kwargs):
        """
        Simple stream/waveform plotter to visualize synthetic waveforms created
        by the solver. Uses ObsPy under the hood to generate a large stream
        and then plots all waveforms together.

        .. note::
            Very simple function to look at waveforms. If you want more
            sophisticated plotting tools, look at Python packages `Pyatoa`
            or `PySEP`

        :type fids: list
        :param fids: list of file ID's to plot
        :type data_format: str
        :param data_format: data format used to determine how to read data files
        :type savefig: str or None
        :param savefig: full path and filename to save the output figure. If
            NoneType, will not save the figure
        """
        from obspy import Stream
        from seisflows.preprocess.default import Default

        # Take advantage of the Default Preprocessing module's read() function
        plotter = Default(data_format=data_format)
        assert(data_format.upper() in plotter._acceptable_data_formats), \
            f"data format must be in {plotter._acceptable_data_formats}"  # NOQA

        st = Stream()
        for fid in fids:
            st += plotter.read(fid)

        st.plot(outfile=savefig, **kwargs)

    def plot2d(self, name=None, parameter=None, cmap=None, savefig=None,
               **kwargs):
        """
        Plot model, gradient or kernels in the PATH.OUTPUT

        :type name: str
        :param name: Name of directory in the output/ directory
        :type parameter: str
        :param parameter: Name of parameter to plot from `name`, e.g., 'vs',
            'vp' etc.
        :type cmap: str
        :param cmap: optional colormap parameter to be passed to Pyplot
        :type savefig: str
        :param savefig: optional name and path of filename to save figure
            to disk
        """
        from seisflows.tools.model import Model

        # Figure out which models/gradients/kernels we can actually plot
        _, output_dir, _ = getpar(key="path_output",
                                  file=self._args.parameter_file,
                                  delim=":")
        # Assuming only models/kernels/gradients have the format *_* in output
        acceptable_names = sorted([
            os.path.basename(_) for _ in glob(os.path.join(output_dir, "*_*"))
        ])
        if name is None:
            print(msg.cli(f"Available models/gradients/kernels",
                          items=sorted(acceptable_names), header="Plot2D")
                  )
            sys.exit(0)
        else:
            assert(name in acceptable_names), (
                f"`seisflows plot 2d` can only plot {acceptable_names}"
            )
        # Grab model_init to use its coordinates
        # name of directory for model_init is defined by solver.specfem.setup()
        base_model = Model(path=os.path.join(output_dir, "MODEL_INIT"))
        assert(base_model.coordinates is not None), \
            f"`MODEL_INIT` does not have any available 2D coordinates"

        # Now read in the actual updated values and update the model
        plot_model = Model(path=os.path.join(output_dir, name))
        plot_model.coordinates = base_model.coordinates
        # plot2d has internal check for acceptable parameter value
        plot_model.plot2d(parameter=parameter, cmap=cmap, show=True,
                          title=f"{name} // {parameter.upper()}", save=savefig)

    def reset(self, choice=None, **kwargs):
        """
        Mid-level function to wrap lower level reset functions

        TODO re-write '_reset_line_search'
        """
        acceptable_args = {"line_search": self._reset_line_search,}

        # Ensure that help message is thrown for empty commands
        if choice not in acceptable_args.keys():
            self._subparser.print_help()
            sys.exit(0)

        acceptable_args[choice](*self._args.args, **kwargs)

    def _inspect_class_that_defined_method(self, name, func, **kwargs):
        """
        Given a function name and generalized module (e.g. solver), inspect
        which of the subclasses actually defined the function. Makes it easier
        to debug/edit source code as the user can quickly determine where
        in the source code they need to look to find the corresponding function.

        https://stackoverflow.com/questions/961048/get-class-that-defined-method

        :type name: str
        :param name: SeisFlows module name
        :type func: str
        :param func: Corresponding method/function name for the given module
        """
        # Dynamically get the correct module and function based on names
        # try:
        #     module = sys.modules[f"seisflows_{name}"]
        # except KeyError:
        #     print(msg.cli(f"SeisFlows has no module: {name}"))
        #     sys.exit(-1)
        # try:
        #     method = getattr(module, func)
        # except AttributeError:
        #     print(msg.cli(f"SeisFlows.{name} has no function: {func}"))
        #     sys.exit(-1)

        parameters = load_yaml(os.path.join(self._args.workdir,
                                            self._args.parameter_file))
        module = custom_import(name, parameters[name])()
        method = getattr(module, func)
        method_name = method.__name__
        if method.__self__:
            classes = [method.__self__.__class__]
        else:
            # Deal with unbound method
            classes = [method.im_class]
        while classes:
            c = classes.pop()
            if method_name in c.__dict__:
                print(f"\n\t{c.__module__}.{c.__name__}.{func}\n")
                return
            else:
                classes = list(c.__bases__) + classes
        print(msg.cli(f"Error matching class for SeisFlows.{name}.{func}"))
        sys.exit(-1)

    def _inspect_module_hierarchy(self, name=None, **kwargs):
        """
        Determine the order of class hierarchy for a given SeisFlows module.

        https://stackoverflow.com/questions/1401661/
                            list-all-base-classes-in-a-hierarchy-of-given-class

        .. rubric::
            $ seisflows print inherit

        :type name: str
        :param name: choice of module, if None, will print hierarchies for all
            modules.
        """
        parameters = load_yaml(os.path.join(self._args.workdir,
                                            self._args.parameter_file))

        items = []
        for NAME in NAMES:
            if name and NAME != name:
                continue
            module = custom_import(NAME, parameters[NAME])()
            item_str = f"{NAME.upper():<12}"
            for i, cls in enumerate(getmro(type(module))[::-1]):
                # The base inheritance is always 'object', skip printing this.
                if i == 0:
                    continue
                item_str += f"> {cls.__name__:<10}"
            items.append(item_str)
        print(msg.cli(items=items, header="seisflows inheritance"))


    def _print_modules(self, package=None, **kwargs):
        """
        Print out available modules in the SeisFlows name space for all
        available packages and modules.

        .. rubric::
            $ seisflows print module

        :type name: str
        :param name: specify an specific module name to list
        :type package: str
        :param package: specify an indivdual package to search
        """
        items = []
        module_dict = return_modules()

        for module_, module_list in module_dict.items():
            if package is not None and module_ != package:
                continue
            items.append(f"- {module_}".expandtabs(tabsize=4))
            for module_ in module_list:
                items.append(f"\t* {module_}".expandtabs(tabsize=4))
        print(msg.cli("'-': module, '*': class", items=items,
                      header="seisflows modules"))

    def _print_tasks(self, **kwargs):
        """
        Simply print out the seisflows.workflow.main() flow variable which
        describes what order workflow functions will be run. Useful for
        filling out the RESUME_FROM and STOP_AFTER parameters.

        .. rubric::
            $ seisflows print flow
        """
        parameters = load_yaml(os.path.join(self._args.workdir,
                                            self._args.parameter_file))
        wf = custom_import("workflow", parameters["workflow"])()
        items = [f"{a+1}: {b.__name__}" for a, b in enumerate(wf.task_list)]
        print(msg.cli(f"Task list for {type(wf)}", items=items,
                      header="seisflows workflow task list"))

    def _print_inheritance(self, name=None, func=None, **kwargs):
        """
        Inspect inheritance hierarchy of classes, methods defined by SeisFlows.
        Useful when developing or debugging, facilitates identification of
        the package top-level.

        USAGE

            seisflows inspect [name] [method]

            To view overall hierarchy for all names in the SeisFlows namespace

                seisflows inspect

            To check the inheritance hierarchy of the 'workflow' module

                seisflows inspect workflow

            To check which class defined a given method, e.g. the 'eval_func'
            method attributed to the solver module

                seisflows inspect solver eval_func

        """
        if func is None:
            self._inspect_module_hierarchy(name, **kwargs)
        else:
            self._inspect_class_that_defined_method(name, func, **kwargs)

    # def _check_model_parameters(self, src=None, **kwargs):
    #     """
    #     Print out the min/max values from one or all of the currently available
    #     models. Useful for checking what models are associated with what part of
    #     the workflow, e.g. evaluate function, evaluate gradient.
    #
    #     :type src: str
    #     :param src: the name of a specific model to check, e.g. 'm_try',
    #         otherwise will check parameters for all models
    #     """
    #     optimize = sys.modules["seisflows_optimize"]
    #     PATH = sys.modules["seisflows_paths"]
    #
    #     avail = glob(os.path.join(PATH.OPTIMIZE, "m_*"))
    #     srcs = [os.path.basename(_) for _ in avail]
    #     if src:
    #         if src not in srcs:
    #             print(msg.cli(f"{src} not in available models: {avail}"))
    #             sys.exit(-1)
    #         srcs = [src]
    #     for tag in srcs:
    #         m = optimize.load(tag)
    #         m.check()
    #
    # def _check_current_iteration(self, **kwargs):
    #     """
    #     Display the current point in the workflow in terms of the iteration
    #     and step count number. Args are not used by allow for a more general
    #     check() function.
    #     """
    #     optimize = sys.modules["seisflows_optimize"]
    #     try:
    #         items = []
    #         ln = optimize.line_search
    #         items.append(f"Iteration: {optimize.iter}")
    #         items.append(f"Step Count: {ln.step_count} / {ln.step_count_max}")
    #         print(msg.cli(items=items))
    #     except AttributeError:
    #         print(msg.cli("OPTIMIZATION module has not been initialized yet, "
    #                       "cannot retrieve iteration or step count values."))
    #         sys.exit(-1)
    #
    # def _check_source_names(self, source_name=None, **kwargs):
    #     """
    #     Sources are tagged by name but also by index in the source names which
    #     can be confusing and usually requires doubling checking. This check
    #     just prints out source names next to their respective index, or if a
    #     source name is requested, provides the index for that
    #
    #     :type source_name: str
    #     :param source_name: name of source to check index, if None will simply
    #         print out all sources
    #     """
    #     try:
    #         source_names = sys.modules["seisflows_solver"].source_names
    #     except FileNotFoundError as e:
    #         print(msg.cli(str(e)))
    #         sys.exit(-1)
    #
    #     if source_name:
    #         print(msg.cli(f"{source_names.index(source_name)}: {source_name}"))
    #     else:
    #         items = []
    #         for i, source_name in enumerate(source_names):
    #             items.append(f"{i:>3}: {source_name}")
    #         print(msg.cli(items=items, header="source names"))
    #
    # def _check_source_index(self, idx=None, **kwargs):
    #     """
    #     Look up source name by index
    #
    #     :type idx: int
    #     :param idx: index of source to look up
    #     """
    #     if idx is None:
    #         self._check_source_names(source_name=None)
    #     else:
    #         solver = sys.modules["seisflows_solver"]
    #         try:
    #             print(msg.cli(f"{idx}: {solver.source_names[int(idx)]}"))
    #         except IndexError:
    #             print(msg.cli(f"idx out of range: {len(solver.source_names)}"))


def return_modules():
    """
    Search for the names of available modules in SeisFlows name space.
    This simple function checks for files with a '.py' extension inside
    each of the sub-directories, ignoring private files like __init__.py.

    :rtype: dict of dict of lists
    :return: a dict with keys matching names and values as dicts for each
        package. nested list contains all the avaialble modules
    """
    module_dict = {}
    for NAME in NAMES:
        module_dict[NAME] = []
        mod_dir = os.path.join(ROOT_DIR, NAME)
        for pyfile in sorted(glob(os.path.join(mod_dir, "*.py"))):
            stripped_pyfile = os.path.basename(pyfile)
            stripped_pyfile = os.path.splitext(stripped_pyfile)[0]
            if not stripped_pyfile.startswith("_"):
                module_dict[NAME].append(stripped_pyfile)

    return module_dict


def main():
    """
    Main entry point into the SeisFlows package is via the SeisFlows class
    """
    sf = SeisFlows()
    sf()


if __name__ == "__main__":
    main()
