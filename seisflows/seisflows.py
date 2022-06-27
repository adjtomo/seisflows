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
"""
import os
import sys
import pickle
import inspect
import logging
import warnings
import argparse
import traceback
import subprocess
from glob import glob
from copy import copy
from IPython import embed

from seisflows.core import Dict, SeisFlowsPathsParameters
from seisflows.config import (config_logger, custom_import, save,
                              NAMES, ROOT_DIR, CFGPATHS)
from seisflows.tools import unix, msg
from seisflows.tools.specfem import (getpar, setpar, getpar_vel_model,
                                     setpar_vel_model)
from seisflows.tools.wrappers import loadyaml


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
                        default=CFGPATHS.PAR_FILE,
                        help=f"Parameters file, default: '{CFGPATHS.PAR_FILE}'")

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
    init = subparser.add_parser(
        "init", help="Initiate working environment",
        description="""Establish a SeisFlows working environment but don't 
        submit the workflow to the system and do not perform variable  error 
        checking. Saves the initial state as pickle files to allow for active 
        environment inspection prior to running 'submit'. Useful for debugging, 
        development and code exploration."""
    )
    # init.add_argument("-c", "--check", action="store_true",
    #                   help="Perform parameter and path checking to ensure that "
    #                        "user-defined parameters are accepatable")
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
    submit.add_argument("-f", "--force", action="store_true",
                        help="Turn off the default parameter precheck")
    submit.add_argument("-s", "--stop_after", default=None, type=str,
                        help="Optional override of the 'STOP_AFTER' parameter")
    # =========================================================================
    resume = subparser.add_parser(
        "resume", help="Re-submit previous workflow to system",
        description="""Resume a previously submitted workflow. Used when 
        an active environment exists in the working directory, and must be 
        submitted to the system again."""
    )
    resume.add_argument("-f", "--force", action="store_true",
                        help="Turn off the default parameter precheck")
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
                         help="Skip the clean and submit precheck statements")
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
    par.add_argument("-r", "--required", action="store_true", default=False,
                     help="Only list parameters which have not been set as a "
                          "default value, typically set with some attention "
                          "catching argument. 'parameter' and 'value' will be "
                          "ignored.")
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
        description="""
Check parameters, state, or values of an active environment

    model     check the min/max values of currently active models tracked by
              optimize. 'seisflows check model [name]' to check specific model.
    iter      Check current interation and step count of workflow
    src       List source names and respective internal indices
    isrc      Check source name for corresponding index
                """,
        help="Check state of an active environment")

    check.add_argument("choice", type=str,  nargs="?",
                       help="Parameter, state, or value to check")
    check.add_argument("args", type=str,  nargs="*",
                       help="Generic arguments passed to check functions")
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
    subparser.add_parser("convert", help="Convert model file format", )
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
    examples.add_argument("run", type=str, nargs="?", default=None,
                          help="Run your choice of example problem")
    examples.add_argument("choice", type=str,  nargs="?", default=None,
                          help="Name of the specific example problem to run")
    # =========================================================================
    # Defines all arguments/functions that expect a sub-argument
    subparser_dict = {"check": check, "par": par, "inspect": inspect,
                      "sempar": sempar, "clean": clean,
                      "restart": restart, "print": print_, "reset": reset,
                      "examples": examples}
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
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        Parse user-defined arguments and establish internal parameters used to
        control which functions execute and how. Instance must be called to
        execute internal functions
        """
        self._parser, self._subparser = sfparser()
        self._paths = None
        self._parameters = None
        self._args = self._parser.parse_args()

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

    def _register_parameters(self, force=True):
        """
        Load the paths and parameters from file into sys.modules, set the
        default parameters if they are missing from the file, and expand all
        paths to absolute pathnames. Also configure the logger.

        .. note::
            This is ideally the FIRST thing that happens everytime SeisFlows
            is initiated. The package cannot do anything without the resulting
            PATH and PARAMETER variables.

        :type force: bool
        :param force: if False, print out a few key parameters and require
            user-input before allowing workflow to be submitted. This is
            usually run before submit and resume, to prevent job submission
            without user evaluation.
        """
        # Check if the filepaths exist
        if not os.path.exists(self._args.parameter_file):
            print(msg.cli(f"SeisFlows parameter file not found: "
                          f"'{self._args.parameter_file}'. Run 'seisflows "
                          f"setup' to create a new parameter file.")
                  )
            sys.exit(-1)

        # Register parameters from the parameter file
        try:
            parameters = loadyaml(self._args.parameter_file)
        except Exception as e:
            print(msg.cli(f"Please check that your parameter file is properly "
                          f"formatted in the YAML format. If you have just run "
                          f"'seisflows configure', you may have some required "
                          f"parameters that will need to be filled out before "
                          f"you can proceed. The error message is:",
                          items=[str(e)], header="parameter file read error",
                          border="="))
            sys.exit(-1)

        # Distribute the paths and parameters internally and separately
        # If we are running seisflows configure, paths will be empty
        try:
            paths = parameters.pop("PATHS")
        except KeyError:
            paths = {}

        # WORKDIR needs to be set here as it's expected by most modules
        if "WORKDIR" not in paths:
            paths["WORKDIR"] = self._args.workdir

        # Parameter file is set here as well so that it can be user-defined
        if "PAR_FILE" not in paths:
            paths["PAR_FILE"] = self._args.parameter_file

        # For submit() and resume(), provide a dialogue to stdout requiring a
        # visual pre-check of parameters before submitting workflow
        # if not force and parameters["PRECHECK"]:
        #     self._precheck_parameters(parameters)

        # Expand all paths to be absolute on the filesystem
        for key, val in paths.items():
            try:
                paths[key] = os.path.expanduser(os.path.abspath(val))
            except TypeError:
                continue

        # Register parameters to sys and internally
        sys.modules["seisflows_parameters"] = Dict(parameters)
        sys.modules["seisflows_paths"] = Dict(paths)
        self._paths = Dict(paths)
        self._parameters = Dict(parameters)

    # def _precheck_parameters(self, parameters):
    #     """
    #     Visually display a list of user-chosen parameters to the User before
    #     proceeding with the _register_parameters command. Allows the User to quickly
    #     determine if workflow parameters have been set correctly
    #
    #     :type parameters: dict
    #     :param parameters: parameters read in from the YAML parameter file
    #     """
    #     items = []
    #     for p in parameters["PRECHECK"]:
    #         try:
    #             items.append(f"{p.upper()}: {parameters[p.upper()]}")
    #         except KeyError:
    #             items.append(f"{p.upper()}: !!! PARAMETER NOT FOUND !!!")
    #     print(msg.cli("Please ensure that the parameters listed below "
    #                   "are set correctly. You can edit this list with "
    #                   "the PRECHECK parameter.", items=items,
    #                   header="seisflows precheck", border="="))
    #     check = input("Continue? (y/[n])\n")
    #     if check != "y":
    #         sys.exit(-1)

    def _register_modules(self, check=True):
        """
        First time setup procedure which loads in the user-chosen modules
        and registers them into sys.modules so that they are globally accessible
        to the program.

        :type check: bool
        :param check: run the check() function for each of the instantiated
            modules, which essentially checks the validity of all the
            user-defined parameters. This is typically wanted, but sometimes
            you don't want to check, e.g., during testing when you know
            some parameters are set incorrectly
        """
        assert(self._paths is not None), (
            f"seisflows._register_parameters() must be run before "
            f"_register_modules()"
        )
        assert(self._parameters is not None), (
            f"seisflows._register_parameters() must be run before "
            f"_register_modules()"
        )
        # Check if current workflow exists on disk, exit so as to not overwrite
        if "OUTPUT" in self._paths and os.path.exists(self._paths.OUTPUT):
            print(msg.cli(
                "Data from previous workflow found in working directory.",
                items=["> seisflows restart: delete data, start new workflow",
                       "> seisflows resume: resume existing workflow"],
                header="warning", border="=")
            )
            sys.exit(-1)

        # Instantiate and register objects
        for name in NAMES:
            sys.modules[f"seisflows_{name}"] = custom_import(name)()

        if check:
            self._check_modules()

        # Bare minimum Module requirements for SeisFlows
        req_modules = ["WORKFLOW", "SYSTEM"]
        for req in req_modules:
            if not hasattr(self._parameters, req):
                print(msg.cli(f"SeisFlows requires modules: {req_modules}."
                              "Please specify these in the parameter file. Use "
                              "'seisflows print module' to determine suitable "
                              "choices.", header="error", border="="))
                sys.exit(-1)

    def _check_modules(self):
        """
        Runs the .check() function on each of the modules, which validates the
        given parameters in a parameter file to ensure that a workflow will not
        break unexpectedly
        """
        errors = []
        for name in NAMES:
            try:
                sys.modules[f"seisflows_{name}"].check()
            except AssertionError as e:
                errors.append(f"{name}: {e}")
        if errors:
            print(msg.cli("seisflows.config module check failed with:",
                          items=errors, header="module check error",
                          border="="))
            sys.exit(-1)

    def _load_modules(self):
        """
        A function to load and check each of the SeisFlows modules,
        re-initiating the SeisFlows environment. All modules are reliant on one
        another so any access to SeisFlows requires loading everything
        simultaneously and in correct order.

        .. note::
            This is similar to config.load() except it doesn't load paths
            and parameters. This allows the User to OVERLOAD the currently
            defined paths and parameters anytime they call 'seisflows resume'
        """
        for NAME in NAMES:
            fid = os.path.join(self._paths.OUTPUT, f"seisflows_{NAME}.p")

            if not os.path.exists(fid):
                print(msg.cli("Not a SeisFlows working directory (no state "
                              "files found). Run 'seisflows init' or "
                              "'seisflows submit' to instantiate a working "
                              "directory.")
                      )
                sys.exit(-1)

            with open(fid, "rb") as f:
                sys.modules[f"seisflows_{NAME}"] = pickle.load(f)
        self._check_modules()

    def setup(self, force=False, **kwargs):
        """
        Initiate a SeisFlows working directory from scratch by establishing a
        template parameter file.

        .. note::
            Future working directory setup functions can be placed here

        :type symlink: bool
        :param symlink: flag to turn on source code symlinking
        :type force: bool
        :param force: flag to force parameter file overwriting
        """
        PAR_FILE = os.path.join(ROOT_DIR, "templates", "parameters.yaml")

        if os.path.exists(self._args.parameter_file):
            if force:
                check = "y"
            else:
                check = input(
                    msg.cli(f"Existing parameter file "
                            f"({self._args.parameter_file}) found. Do you "
                            f"wish to overwrite with a blank file? (y/[n])"
                            ))
            if check == "y":
                unix.rm(self._args.parameter_file)
            else:
                sys.exit(0)

        unix.cp(PAR_FILE, self._args.workdir)
        print(msg.cli(f"creating parameter file: {self._args.parameter_file}"))

    def configure(self, absolute_paths=False, **kwargs):
        """
        Dynamically generate the parameter file by writing out docstrings and
        default values for each of the SeisFlows module parameters.
        This function writes files manually, consistent with the .yaml format.

        :type absolute_paths: bool
        :param absolute_paths: if True, expand pathnames to absolute paths,
            else if False, use path names relative to the working directory.
            Defaults to False, uses relative paths.
        """
        print(msg.cli(f"filling {self._args.parameter_file} w/ default values"))
        self._register_parameters(force=True)

        # Check if the User set turn off any modules (if None, dont instantiate)
        names = copy(NAMES)
        for name, choice in self._parameters.items():
            if choice is None:
                names.remove(name.lower())

        # Need to attempt importing all modules before we access their par/paths
        for NAME in NAMES:
            sys.modules[f"seisflows_{NAME}"] = custom_import(NAME)()

        # System defines foundational directory structure required by other
        # modules. Don't validate the parameters because they aren't yet set
        sys.modules["seisflows_system"].required.validate(paths=True,
                                                          parameters=False)

        # If writing to parameter file fails for any reason, the file will be
        # mangled, create a temporary copy that can be re-instated upon failure
        temp_par_file = f".{self._args.parameter_file}"
        unix.cp(self._args.parameter_file, temp_par_file)

        try:
            # Paths are collected for each but written at the end
            seisflows_paths = {}
            with open(self._args.parameter_file, "a") as f:
                for name in names:
                    req = sys.modules[f"seisflows_{name}"].required
                    seisflows_paths.update(req.paths)

                    # Write the docstring header and then the parameters in YAML
                    msg.write_par_file_header(f, req.parameters, name)
                    msg.write_par_file_paths_pars(f, req.parameters)

                # Write the paths in the same format as parameters
                msg.write_par_file_header(f, seisflows_paths, name="PATHS")
                f.write("PATHS:\n")

                # If requested, set the paths relative to the current dir
                if not absolute_paths:
                    for key, attrs in seisflows_paths.items():
                        if attrs["default"] is not None:
                            seisflows_paths[key]["default"] = os.path.relpath(
                                                               attrs["default"])
                msg.write_par_file_paths_pars(f, seisflows_paths, indent=4)
        # General error catch as anything can happen here
        except Exception as e:
            unix.rm(self._args.parameter_file)
            unix.cp(temp_par_file, self._args.parameter_file)
            print(msg.cli(text="seisflows configure traceback", header="error"))
            print(traceback.format_exc())
            sys.exit(-1)
        else:
            unix.rm(temp_par_file)

    def init(self, **kwargs):
        """
        Establish a SeisFlows working environment on disk. Instantiates a
        working state in memory (sys.modules) and then writes this state as \
        pickle files to the OUTPUT directory for User inspection and debug
        purposes.
        """
        unix.mkdir(self._args.workdir)
        unix.cd(self._args.workdir)

        self._register_parameters(force=True)
        self._register_modules(check=True)

        save(path=self._paths.OUTPUT)

        # Ensure that all parameters and paths that need to be instantiated
        # are present in sys modules
        for NAME in NAMES:
            sys.modules[f"seisflows_{NAME}"].required.validate()

        print(msg.cli(f"instantiating SeisFlows working state in directory: "
                      f"{CFGPATHS.OUTPUTDIR}"))

    def submit(self, stop_after=None, force=False, **kwargs):
        """
        Main SeisFlows execution command. Submit the SeisFlows workflow to
        the chosen system, and execute seisflows.workflow.main(). Will create
        the working directory and any required paths and ensure that all
        required paths exist.

        :type stop_after: str
        :param stop_after: allow the function to overwrite the 'STOP_AFTER'
            parameter in the parameter file, which dictates how far the workflow
            will proceed until stopping. Must match flow function names in
            workflow.main()
        :type force: bool
        :param force: if True, turns off the parameter precheck and
            simply submits the workflow
        """
        unix.mkdir(self._args.workdir)
        unix.cd(self._args.workdir)

        # Ensure that the 'RESUME_FROM' parameter is not set, in case of restart
        self.par(parameter="resume_from", value="", skip_print=True)
        if stop_after is not None:
            self.par(parameter="stop_after", value=stop_after, skip_print=True)

        # Read in the Parameter file and set parameters into sys.modules.
        self._register_parameters(force=force)
        # self._check_required_paths()
        self._register_modules(check=True)

        # Set logger to print to stdout and write to a file
        config_logger(level=self._parameters.LOG_LEVEL,
                      verbose=self._parameters.VERBOSE,
                      filename=self._paths.LOGFILE)

        # Submit workflow.main() to the system
        system = sys.modules["seisflows_system"]
        system.submit()

    # def _check_required_paths(self):
    #     """
    #     If the User provides certain paths to the program, they MUST exist.
    #     This function simply checks these required paths and throws a sys exit
    #     if any of them does not exist
    #     """
    #     # A list of paths that need to exist if provided by user
    #     REQ_PATHS = ["SPECFEM_BIN", "SPECFEM_DATA", "MODEL_INIT", "MODEL_TRUE",
    #                  "DATA", "LOCAL", "MASK"]
    #
    #     # Check that all required paths exist before submitting workflow
    #     paths_dont_exist = []
    #     for key in REQ_PATHS:
    #         if key in self._paths:
    #             # If a required path is given (not None) and doesnt exist, exit
    #             if self._paths[key] and not os.path.exists(self._paths[key]):
    #                 paths_dont_exist.append(f"{key}: {self._paths[key]}")
    #     if paths_dont_exist:
    #         print(msg.cli("The following paths do not exist but need to:",
    #                       items=paths_dont_exist, header="path error",
    #                       border="="))
    #         sys.exit(-1)

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
            # CFGPATHS defines the outermost directory structure of SeisFlows
            # We safeguard below against deleting the parameter file
            items = []
            for fid_ in CFGPATHS.values():
                for fid in glob(os.path.join(self._args.workdir, fid_)):
                    # Safeguards against deleting files that should not be dltd
                    try:
                        assert("yaml" not in fid)
                        assert(not os.path.islink(fid))
                        unix.rm(fid)
                        items.append(f"- deleting file/folder: {fid}")
                    except AssertionError:
                        items.append(f"+ skipping over: {fid}")
                        continue
            print(msg.cli(items=items, header="clean", border="="))

    def resume(self, stop_after=None, resume_from=None, force=False,
               **kwargs):
        """
        Resume a previously started workflow by loading the module pickle files
        and submitting the workflow from where it left off.

        :type stop_after: str
        :param stop_after: allow the function to overwrite the 'STOP_AFTER'
            parameter in the parameter file, which dictates how far the workflow
            will proceed until stopping. Must match flow function names in
            workflow.main()
        :type resume_from: str
        :param resume_from: allow the function to overwrite the 'RESUME_FROM'
            parameter in the parameter file, which dictates which function the
            workflow starts from, must match the flow functions given in
            workflow.main()
        :type force: bool
        :param force: if True, turns off the parameter precheck and
            simply submits the workflow
        """
        if stop_after is not None:
            self.par(parameter="STOP_AFTER", value=stop_after, skip_print=True)
        if resume_from is not None:
            self.par(parameter="RESUME_FROM", value=resume_from,
                     skip_print=True)

        self._register_parameters(force=force)
        self._load_modules()

        # Set logger to print to stdout and write to a file
        config_logger(level=self._parameters.LOG_LEVEL,
                      verbose=self._parameters.VERBOSE,
                      filename=self._paths.LOGFILE)

        system = sys.modules["seisflows_system"]
        system.submit()

    def restart(self, force=False, **kwargs):
        """
        Restart simply means clean the workding dir and submit a new workflow.

        :type force: bool
        :param force: ignore the warning check that precedes the clean() 
            function, useful if you don't want any input messages popping up
        """
        self.clean(force=force)
        self.submit(force=force)

    def debug(self, **kwargs):
        """
        Initiate an IPython debugging environment to explore the currently
        active SeisFlows environment. Reloads the system modules in an
        interactive environment allowing exploration of the package space.
        Does not allow stepping through of code (not a breakpoint).
        """
        self._register_parameters(force=True)
        self._load_modules()

        # Distribute modules to common names for easy access during debug mode
        PATH = sys.modules["seisflows_paths"]
        PAR = sys.modules["seisflows_parameters"]
        system = sys.modules["seisflows_system"]
        preprocess = sys.modules["seisflows_preprocess"]
        solver = sys.modules["seisflows_solver"]
        postprocess = sys.modules["seisflows_postprocess"]
        optimize = sys.modules["seisflows_optimize"]
        workflow = sys.modules["seisflows_workflow"]

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
                sys.exit(-1)

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

    def par(self, parameter, value=None, skip_print=False, required=False,
            **kwargs):
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
        :type required: bool
        :param required: Only list parameters which have not been set as a
            default value, 'parameter' and 'value' will be ignored.
        """
        if not os.path.exists(self._args.parameter_file):
            sys.exit(f"\n\tparameter file '{self._args.parameter_file}' "
                     f"does not exist\n")

        if parameter is None and not required:
            self._subparser.print_help()
            sys.exit(0)
        elif required:
            self._par_required()
            sys.exit(0)

        # SeisFlows parameter file dictates upper-case parameters
        parameter = parameter.upper()
        if isinstance(value, str) and value.lower() == "none":
            warnings.warn("to set values to NoneType, use 'null' not 'none'",
                          UserWarning)

        # Use the specfem tool to grab related information
        try:
            key, cur_val, i = getpar(key=parameter,
                                     file=self._args.parameter_file,
                                     delim=":")
        except KeyError:
            print(msg.cli(f"'{parameter}' not found in "
                          f"{self._args.parameter_file}"))
            sys.exit(-1)

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

    def _par_required(self):
        """
        Only list parameters which have not been set as a default value.
        Filled in with default values defined in SeisFlowsPathParameters
        """
        sf = SeisFlowsPathsParameters
        with open(self._args.parameter_file, "r") as f:
            lines = f.readlines()
            for check in [sf.default_par, sf.default_path]:
                print(f"{check}\n{'='*len(check)}")
                for line in lines:
                    if check in line:
                        print(f"\t{line.split(':')[0].strip()}")

    def examples(self, run=None, choice=None, **kwargs):
        """
        List or run a SeisFlows example problem

        USAGE

            seisflows examples [run] [choice]

            To list available examples:

                seisflows examples

            To run a specific example (this is the same as 'python example.py')

                seisflows examples run 1

        :type run: bool
        :param run: if True, run an example of choice `choice`
        :type choice: str
        :param choice: The choice of example, must match the given tag or file
            name that is assigned to it
        """
        examples_dir = os.path.join(ROOT_DIR, "scripts", "examples")
        examples_list = []
        example_names = sorted(glob(os.path.join(examples_dir, "ex*.py")))

        for i, fid in enumerate(example_names):
            example_name = os.path.splitext(os.path.basename(fid))[0]
            examples_list.append((i+1, example_name, fid))

        arg1, arg2 = None, None
        if run:
            # Case 1: seisflows examples 1 OR seisflows examples ex1_...
            if choice is None:
                arg1 = run
                arg2 = ""
            # Case 2: seisflows examples run 1 OR seisflows examples run ex1_...
            elif run in ["run", "setup"]:
                arg1 = choice
                arg2 = f" {run}"  # space so that we do $ python ex.py run
        if arg1:
            # Allow for matching against index (int) and name (str)
            try:
                arg1 = int(arg1)
            except ValueError:
                pass

            for ex_tup in examples_list:
                j, exname, fid = ex_tup
                if arg1 in [j, exname]:
                    print(f"{run.capitalize()} example: {exname}")
                    subprocess.run(f"python {fid}{arg2}", shell=True,
                                   check=False)
                    return

        # Default behavior is to just print this help dialogue
        items = [f"{j}: {exname}" for j, exname, fid in examples_list]
        print(msg.cli("Example options where <name_or_idx> is either the "
                      "example name or corresponding index, provided below.",
                      items=[
            "'seisflows examples <name_or_idx>': print example description",
            "'seisflows examples setup <name_or_idx>': setup example but "
            "don't run workflow",
            "'seisflows examples run <name_or_idx>': setup and run example"
        ],
            header="seisflows examples"
        ))
        print(msg.cli(items=items))

    def check(self, choice=None, **kwargs):
        """
        Check parameters, state or values  of an active SeisFlows environment.
        Type 'seisflows check --help' for a detailed help message.

        :type choice: str
        :param choice: underlying sub-function to choose
        """
        acceptable_args = {"model": self._check_model_parameters,
                           "iter": self._check_current_iteration,
                           "src": self._check_source_names,
                           "isrc": self._check_source_index}

        # Ensure that help message is thrown for empty commands
        if choice not in acceptable_args.keys():
            self._subparser.print_help()
            sys.exit(0)

        self._register_parameters(force=True)
        self._load_modules()
        acceptable_args[choice](*self._args.args, **kwargs)

    def print(self, choice=None, **kwargs):
        """
        Print information relating to an active SeisFlows environment.
        Type 'seisflows check --help' for a detailed help message.

        :type choice: str
        :param choice: underlying sub-function to choose
        """
        acceptable_args = {"modules": self._print_modules,
                           "flow": self._print_flow,
                           "inherit": self._print_inheritance}

        # Ensure that help message is thrown for empty commands
        if choice not in acceptable_args.keys():
            self._subparser.print_help()
            sys.exit(0)

        acceptable_args[choice](*self._args.args, **kwargs)

    def reset(self, choice=None, **kwargs):
        """
        Mid-level function to wrap lower level reset functions
        """
        acceptable_args = {"line_search": self._reset_line_search,}

        # Ensure that help message is thrown for empty commands
        if choice not in acceptable_args.keys():
            self._subparser.print_help()
            sys.exit(0)

        self._register_parameters(force=True)
        self._load_modules()
        acceptable_args[choice](*self._args.args, **kwargs)

    def convert(self, name, path=None, **kwargs):
        """
        Convert a model in the OUTPUT directory between vector to binary
        representation. Kwargs are passed through to solver.save()

        USAGE

            seisflows convert [name] [path] [**kwargs]

            To convert the vector model 'm_try' to binary representation in the
            output directory

                seisflows convert m_try

        :type name: str
        :param name: name of the model to convert, e.g. 'm_try'
        :type path: str
        :param path: path and file id to save the output model. if None, will
            default to saving in the output directory under the name of the
            model
        """
        self._load_modules(force=True)

        solver = sys.modules["seisflows_solver"]
        optimize = sys.modules["seisflows_optimize"]
        PATH = sys.modules["seisflows_paths"]

        if path is None:
            path = os.path.join(PATH.OUTPUT, name)
        if os.path.exists(path):
            print(msg.cli("The following file exists and will be overwritten. "
                          "Please rename or move this file and re-try: {path}"))
            sys.exit(-1)

        solver.save(solver.split(optimize.load(name)), path=path, **kwargs )

    @staticmethod
    def _inspect_class_that_defined_method(name, func, **kwargs):
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
        try:
            module = sys.modules[f"seisflows_{name}"]
        except KeyError:
            print(msg.cli(f"SeisFlows has no module: {name}"))
            sys.exit(-1)
        try:
            method = getattr(module, func)
        except AttributeError:
            print(msg.cli(f"SeisFlows.{name} has no function: {func}"))
            sys.exit(-1)

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

    @staticmethod
    def _inspect_module_hierarchy(name=None, **kwargs):
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
        items = []
        for NAME in NAMES:
            if name and NAME != name:
                continue
            module = sys.modules[f"seisflows_{NAME}"]
            item_str = f"{NAME.upper():<12}"
            for i, cls in enumerate(inspect.getmro(type(module))[::-1]):
                item_str += f"> {cls.__name__:<10}"
            items.append(item_str)
        print(msg.cli(items=items, header="seisflows inheritance"))

    def _reset_line_search(self, **kwargs):
        """
        Reset the machinery of the line search. This is useful for if a line
        search fails or stagnates but the User does not want to re-run the
        entire iteration. They can reset the line search and resume the workflow
        from the line search step

        The following rubric details how you might use this from command line:

        .. rubric::
            $ seisflows reset line_search
            $ seisflows par resume_from line_search
            $ seisflows resume_from -f
        """
        optimize = sys.modules["seisflows_optimize"]
        workflow = sys.modules["seisflows_workflow"]
        
        current_step = optimize.line_search.step_count
        optimize.line_search.reset()

        # Manually set step count back to 0, this usually happens in
        # optimize.finalize_search()
        optimize.line_search.step_count = 0

        print(msg.cli(f"resetting line search machinery. step count: "
                      f"{current_step} -> {optimize.line_search.step_count }"))
        workflow.checkpoint()

    def _print_modules(self, name=None, package=None, **kwargs):
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
        print(msg.cli("'+': package, '-': module, '*': class", items=items,
                      header="seisflows modules"))

    def _print_flow(self, **kwargs):
        """
        Simply print out the seisflows.workflow.main() flow variable which
        describes what order workflow functions will be run. Useful for
        filling out the RESUME_FROM and STOP_AFTER parameters.

        .. rubric::
            $ seisflows print flow
        """
        self._register_parameters(force=True)
        self._load_modules()

        workflow = custom_import("workflow")()
        flow = workflow.main(return_flow=True)
        items = [f"{a+1}: {b.__name__}" for a, b in enumerate(flow)]
        print(msg.cli(f"Flow arguments for {type(workflow)}", items=items,
                      header="seisflows workflow main"))

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
        self._register_parameters(force=True)
        self._load_modules()
        if func is None:
            self._inspect_module_hierarchy(name, **kwargs)
        else:
            self._inspect_class_that_defined_method(name, func, **kwargs)

    def _check_model_parameters(self, src=None, **kwargs):
        """
        Print out the min/max values from one or all of the currently available
        models. Useful for checking what models are associated with what part of
        the workflow, e.g. evaluate function, evaluate gradient.

        :type src: str
        :param src: the name of a specific model to check, e.g. 'm_try', 
            otherwise will check parameters for all models
        """
        optimize = sys.modules["seisflows_optimize"]
        PATH = sys.modules["seisflows_paths"]

        avail = glob(os.path.join(PATH.OPTIMIZE, "m_*"))
        srcs = [os.path.basename(_) for _ in avail]
        if src:
            if src not in srcs:
                print(msg.cli(f"{src} not in available models: {avail}"))
                sys.exit(-1)
            srcs = [src]
        for tag in srcs:
            m = optimize.load(tag)
            optimize.check_model(m, tag)

    def _check_current_iteration(self, **kwargs):
        """
        Display the current point in the workflow in terms of the iteration
        and step count number. Args are not used by allow for a more general
        check() function.
        """
        optimize = sys.modules["seisflows_optimize"]
        try:
            items = []
            ln = optimize.line_search
            items.append(f"Iteration: {optimize.iter}")
            items.append(f"Step Count: {ln.step_count} / {ln.step_count_max}")
            print(msg.cli(items=items))
        except AttributeError:
            print(msg.cli("OPTIMIZATION module has not been initialized yet, "
                          "cannot retrieve iteration or step count values."))
            sys.exit(-1)

    def _check_source_names(self, source_name=None, **kwargs):
        """
        Sources are tagged by name but also by index in the source names which
        can be confusing and usually requires doubling checking. This check
        just prints out source names next to their respective index, or if a
        source name is requested, provides the index for that

        :type source_name: str
        :param source_name: name of source to check index, if None will simply
            print out all sources
        """
        try:
            source_names = sys.modules["seisflows_solver"].source_names
        except FileNotFoundError as e:
            print(msg.cli(str(e)))
            sys.exit(-1)

        if source_name:
            print(msg.cli(f"{source_names.index(source_name)}: {source_name}"))
        else:
            items = []
            for i, source_name in enumerate(source_names):
                items.append(f"{i:>3}: {source_name}")
            print(msg.cli(items=items, header="source names"))

    def _check_source_index(self, idx=None, **kwargs):
        """
        Look up source name by index

        :type idx: int
        :param idx: index of source to look up
        """
        if idx is None:
            self._check_source_names(source_name=None)
        else:
            solver = sys.modules["seisflows_solver"]
            try:
                print(msg.cli(f"{idx}: {solver.source_names[int(idx)]}"))
            except IndexError:
                print(msg.cli(f"idx out of range: {len(solver.source_names)}"))


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
