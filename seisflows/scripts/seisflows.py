#!/usr/bin/env python
"""
The main entry point to the SeisFlows package. A high-level command line tool 
that facilitates interface with the underlying SeisFlows environment, package.
"""
import os
import sys
import inspect
import warnings
import argparse
import subprocess
from glob import glob
from textwrap import wrap
from seisflows.tools import unix, tools
from seisflows.tools.tools import loadyaml, loadpy, parse_null
from seisflows.config import (init_seisflows, tilde_expand, Dict, custom_import,
                              NAMES, PACKAGES, ROOT_DIR)


def main_parser():
    """
    Get User defined arguments, or assign defaults

    :rtype: argparse.ArgumentParser()
    :return: User defined or default arguments
    """
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
These are common SeisFlows commands and their intended usages:

initiate a working directory from scratch:
    setup       Initiate a blank parameter file, symlink the source code
    modules     List out available modules names based on file availability
    configure   Set parameter file with docstrings and defaults based on modules
    init        Initiate a SeisFlows working environment.

submit a SeisFlows workflow to the system:
    submit      Main execution command, submit initial workflow to system
    resume      Re-submit a previously submitted workflow
    clean       Clean up working directory, leave parameter file
    
view / manipulate an active environment (see also: seisflows [command] --help):
    par         View and edit values in the parameter file
    check       Check parameters, state, or values of an active environment
    convert     Convert a model from binary to numpy array or vice-versa
    reset       Reset chosen underlying machinery
    
debug and development tools:
    inspect     Inspect inheritance hierarchy and method ownership
    debug       Start an IPython debugger to explore an active working state
    edit        Edit a given source code file with a preferred text editor
    help        Print out docstrings for methods listed here
""",
    )

    # Positional arguments
    parser.add_argument("command", type=str, nargs="?",
                        help="A task for SeisFlows to perform")
    parser.add_argument("arguments", type=str, nargs="*",
                        help="Optional arguments for given command")

    # Optional parameters
    parser.add_argument("-w", "--workdir", nargs="?", default=os.getcwd(),
                        help="The SeisFlows working directory, default: cwd")
    parser.add_argument("-p", "--parameter_file", nargs="?",
                        default="parameters.yaml",
                        help="Parameters file, default: 'parameters.yaml'")
    parser.add_argument("--path_file", nargs="?", default="paths.py",
                        help="Legacy path file, default: 'paths.py'")


    return parser


class SeisFlows:
    """
    The main entry point to the SeisFlows package. Responsible for setting up
    or re-creating a SeisFlows enviornment, (re-)submitting workflows, or
    inspecting, manipulating or viewing a live environment via command line
    arguments.

    .. note::
        Almost every modules requires loading of other modules, i.e. to run
        any checks we must load the entire SeisFlows environment, which is slow
        but provides the most flexibility when accessing internal information
    """
    def __init__(self):
        """
        Parse user-defined arguments, call internal method with given arguments
        """
        self._parser = main_parser()
        self._args = self._parser.parse_args()

        # To be set by _register()
        self._paths = None
        self._parameters = None

        if self._args.command is None:
            self._parser.print_help()
            sys.exit(0)

        if not hasattr(self, self._args.command):
            sys.exit(f"\n\tno matching function '{self._args.command}'"
                     f"\n\ttype 'seisflows' for available functions\n")
        try:
            getattr(self, self._args.command)(*self._args.arguments)
        except TypeError as e:
            # Missing required parameters will throw TypeErrors, ignore the
            # traceback information and simply throw the formatted error message
            sys.exit(f"\n\tseisflows {e}\n")

    @property
    def _public_methods(self):
        """
        Return a list of all public methods within this class.

        .. warning::
            Only methods that can be called via the command line should be
            public, all other methods and attributes should be private.
        """
        return [_ for _ in dir(self) if not _.startswith("_")]

    def _register(self, precheck=True):
        """
        Load the paths and parameters from file into sys modules, set the
        default parameters if they are missing from the file, and expand all
        paths to absolute pathnames.

        :type precheck: bool
        :param precheck: print out a few key parameters and require user-input
            before allowing workflow to be submitted. This is usually run before
            submit and resume, to prevent job submission without evaluation.
        """
        # Check if the filepaths exist
        if not os.path.exists(self._args.parameter_file):
            sys.exit(f"\n\tSeisFlows parameter file not found: "
                     f"{self._args.parameter_file}\n")

        # Register parameters from the parameter file
        if self._args.parameter_file.endswith(".yaml"):
            parameters = loadyaml(self._args.parameter_file)
            try:
                paths = parameters["PATHS"]
                parameters.pop("PATHS")
            except KeyError:
                paths = {}
        #  Allow for legacy .py parameter file naming
        elif self._args.parameter_file.endwith(".py"):
            warnings.warn(".py parameter and path files are deprecated in "
                          "favor of a .yaml parameter file. Please consider "
                          "switching as the use of legacy .py files may have "
                          "unintended consequences at runtime",
                          DeprecationWarning)

            if not os.path.exists(self._args.paths_file):
                sys.exit(f"\n\tLegacy parameter file requires corresponding "
                         f"path file\n")
            parameters = loadpy(self._args.parameter_file)
            paths = loadpy(self._args.path_file)
        else:
            raise TypeError(f"Unknown file format for "
                            f"{self._args.parameter_file}, file must be "
                            f"'.yaml' (preferred) or '.py' (legacy)")

        # WORKDIR needs to be set here as it's expected by most modules
        if "WORKDIR" not in paths:
            paths["WORKDIR"] = self._args.workdir

        # For submit() and resume(), provide a dialogue to stdout requiring a
        # visual pre-check of parameters before submitting workflow
        if precheck and parameters["PRECHECK"]:
            print("\n\tSEISFLOWS PARAMETER CHECK"
                  "\n\t=========================\n")
            for par in parameters["PRECHECK"]:
                par = par.upper()
                try:
                    print(f"\t{par}: {parameters[par]}")
                except KeyError:
                    print(f"\t{par}: !!! PARAMETER NOT FOUND !!!")
            print("\n")
            check = input("\tContinue? (y/[n]): ")
            if check != "y":
                sys.exit(-1)

        # Register parameters to sys, ensure they meet standards of the package
        # parameters = parse_null(parameters)
        sys.modules["seisflows_parameters"] = Dict(parameters)

        # Register paths to sys, expand to relative paths to absolute, drop null
        paths = tilde_expand(parse_null(paths))
        paths = {key: os.path.abspath(path) for key, path in paths.items()}
        sys.modules["seisflows_paths"] = Dict(paths)

        self._paths = paths
        self._parameters = parameters

    def _load_modules(self, **kwargs):
        """
        A function to load and check each of the SeisFlows modules,
        re-initiating the SeisFlows environment. All modules are reliant on one
        another so any access to SeisFlows requires loading everything
        simultaneously.
        """
        self._register(**kwargs)

        # Working directory should already have been created by submit()
        unix.cd(self._args.workdir)

        # Reload objects from Pickle files
        for name in NAMES:
            fullfile = os.path.join(self._args.workdir, "output",
                                    f"seisflows_{name}.p")

            if not os.path.exists(fullfile):
                sys.exit(f"\n\tNot a SeisFlows working directory, state file "
                         f"not found:\n\t{fullfile}\n")

            sys.modules[f"seisflows_{name}"] = tools.loadobj(fullfile)

        # Check parameters so that default values are present
        for name in NAMES:
            sys.modules[f"seisflows_{name}"].check()

    def help(self, choice=None):
        """
        Help messages regarding available SeisFlows tasks. Prints out the
        docstrings for each of the public methods available in this class.

        :type choice: str
        :param choice: if not None, will request an indvidual help message
        """
        if choice:
            try:
                docstr = getattr(self, choice).__doc__
                print(f"\n> seisflows {choice}\n"
                      f"{docstr}")
            except AttributeError:
                sys.exit(f"\n\tseisflows help has no entry '{choice}'\n")
        else:
            for meth in self._public_methods:
                print(f"\n> seisflows {meth}")
                print(getattr(self, meth).__doc__)

    def modules(self):
        """
        Search for the names of available modules in SeisFlows name space.
        This simple function checks for files with a '.py' extension inside
        each of the sub-directories, ignoring private files like __init__.py.
        """
        REPO_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

        for name in NAMES:
            print(f"\n\t{name.upper()}")
            for package in PACKAGES:
                print(f"\t\t{package}")
                mod_dir = os.path.join(REPO_DIR, package, name)
                for pyfile in sorted(glob(os.path.join(mod_dir, "*.py"))):
                    stripped_pyfile = os.path.basename(pyfile)
                    stripped_pyfile = os.path.splitext(stripped_pyfile)[0]
                    if not stripped_pyfile.startswith("_"):
                        print(f"\t\t\t{os.path.basename(stripped_pyfile)}")
        print("\n")

    def setup(self):
        """
        Initiate a SeisFlows working directory from scratch; establish a
        template parameter file and symlink the source code for easy access
        """
        PAR_FILE =  os.path.join(ROOT_DIR, "templates", "parameters.yaml")
        REPO_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

        if os.path.exists(self._args.parameter_file):
            print(f"\n\tParameter file '{self._args.parameter_file}' "
                  f"already exists\n")
            check = input(f"\tOverwrite with blank file? (y/[n]): ")
            if check == "y":
                unix.rm(self._args.parameter_file)
        unix.cp(PAR_FILE, self._args.workdir)

        # Symlink the source code for easy access to repo
        src_code = os.path.join(self._args.workdir, "source_code")
        if not os.path.exists(src_code):
            unix.mkdir(src_code)
            for package in PACKAGES:
                unix.ln(os.path.join(REPO_DIR, package), src_code)

    def configure(self):
        """
        Dynamically generate the parameter file by writing out docstrings and
        default values for each of the SeisFlows module parameters.
        This function writes files manually, consistent with the .yaml format.
        """
        self._register(precheck=False)

        # Set some re-usable strings to provide a consistent look
        BLANK = "\n"
        COMMENT = "#\n"
        DIVIDER = f"# {'=' * 78}\n"
        HEADER = "# {}\n"
        UNDERLINE = f"# {'-' * 25}\n"
        TAB = "    "

        HEADER_TOP = BLANK + DIVIDER + COMMENT + HEADER + UNDERLINE + COMMENT
        HEADER_BOT = COMMENT + DIVIDER

        # Establish the paths and parameters provided by the user
        if not self._args.parameter_file.endswith(".yaml"):
            sys.exit(f"\n\tseisflows configure only applicable to .yaml "
                     f"parameter files\n")

        # Need to attempt importing all modules before we access any of them
        for name in NAMES:
            sys.modules[f"seisflows_{name}"] = custom_import(name)()

        # System defines foundational directory structure required by other
        # modules. Don't validate the parameters because they aren't yet set
        sys.modules["seisflows_system"].required.validate(paths=True,
                                                          parameters=False)

        # Paths are collected from each module but written at the end together
        seisflows_paths = {}

        # If writing to parameter file fails for any reason, the file will be
        # mangled, create a temporary copy that can be re-instated upon failure
        temp_par_file = f".{self._args.parameter_file}"
        unix.cp(self._args.parameter_file, temp_par_file)
        try:
            with open(self._args.parameter_file, "a") as f:
                for name in NAMES:
                    req = sys.modules[f"seisflows_{name}"].required
                    seisflows_paths.update(req.paths)

                    # Write the docstring header with all parameters
                    f.write(HEADER_TOP.format(name.upper()))
                    for key, attrs in req.parameters.items():
                        f.write(f"# {key} ({attrs['type']}):\n")
                        # Ensure that header lines are no more than 80 char
                        docstrs = wrap(attrs["docstr"], width=77 - len(TAB),
                                       break_long_words=False)
                        for line, docstr in enumerate(docstrs):
                            f.write(f"#{TAB}{docstr}\n")
                    f.write(HEADER_BOT)

                    # Write parameters in a YAML consistent format
                    for key, attrs in req.parameters.items():
                        # Lists need to be treated differently
                        if isinstance(attrs["default"], list):
                            f.write(f"{key}:\n")
                            for val in attrs["default"]:
                                f.write(f"{TAB}- {val}\n")
                        else:
                            # Yaml writes None as 'null'
                            if attrs["default"] is None:
                                attrs["default"] = "null"
                            f.write(f"{key}: {attrs['default']}\n")

                # Write the paths in the same format as parameters
                f.write(HEADER_TOP.format("PATHS"))
                for key, attrs in seisflows_paths.items():
                    # Ensure that header lines are no more than 80 char
                    docstr_ = f"{key}: {attrs['docstr']}"
                    docstrs = wrap(docstr_, width=77 - len(TAB),
                                   break_long_words=False)
                    for line, docstr in enumerate(docstrs):
                        f.write(f"# {docstr}\n")
                f.write(HEADER_BOT)

                f.write("PATHS:\n")
                for key, attrs in seisflows_paths.items():
                    if attrs["default"] is None:
                        default = "null"
                    else:
                        default = attrs["default"]
                    f.write(f"{TAB}{key}: {default}\n")
        except Exception as e:
            # General error catch as anything can happen here
            unix.rm(self._args.parameter_file)
            unix.cp(temp_par_file, self._args.parameter_file)
            sys.exit(f"\n\tseisflows configure failed with exception:\n\t{e}\n")
        else:
            unix.rm(temp_par_file)

    def init(self):
        """
        Establish a SeisFlows working environment without error checking.
        Save the initial state as pickle files for environment inspection.
        Useful for debugging, development and code exploration purposes.
        """
        self._register(precheck=False)

        unix.mkdir(self._args.workdir)
        unix.cd(self._args.workdir)

        # Submit workflow
        init_seisflows()
        workflow = sys.modules["seisflows_workflow"]
        workflow.checkpoint()

    def submit(self):
        """
        Main SeisFlows execution command. Submit the SeisFlows workflow to
        the chosen system, and execute seisflows.workflow.main(). Will create
        the working directory and any required paths and ensure that all
        required paths exist.
        """
        self._register(precheck=True)

        # A list of paths that need to exist if provided by user
        REQ_PATHS = ["SPECFEM_BIN", "SPECFEM_DATA", "MODEL_INIT", "MODEL_TRUE",
                     "DATA", "LOCAL", "MASK"]

        # Check that all required paths exist before submitting workflow
        paths_dont_exist = []
        for key in REQ_PATHS:
            if key in self._paths and not os.path.exists(self._paths[key]):
                paths_dont_exist.append(self._paths[key])
        if paths_dont_exist:
            print("\nThe following paths do not exist:\n")
            for path_ in paths_dont_exist:
                print(f"\t{path_}")
            print("\n")
            sys.exit()

        unix.mkdir(self._args.workdir)
        unix.cd(self._args.workdir)

        # Submit workflow.main() to the system
        init_seisflows()
        workflow = sys.modules["seisflows_workflow"]
        system = sys.modules["seisflows_system"]
        system.submit(workflow)

    def clean(self):
        """
        Clean the SeisFlows working directory except for the parameter file.
        """
        check = input("\n\tThis will remove all workflow objects, leaving only "
                      "the parameter file.\n\tAre you sure you want to clean? "
                      "(y/[n]): ")

        if check == "y":
            for fid in glob(os.path.join(self._args.workdir, "output*")):
                unix.rm(fid)
            for fid in glob(os.path.join(self._args.workdir, "*log*")):
                unix.rm(fid)
            unix.rm(os.path.join(self._args.workdir, "scratch"))

    def resume(self):
        """
        Resume a previously started workflow by loading the module pickle files
        and submitting the workflow from where it left off.
        """
        self._load_modules(precheck=True)

        workflow = sys.modules["seisflows_workflow"]
        system = sys.modules["seisflows_system"]

        system.submit(workflow)

    def restart(self):
        """
        Restart simply means clean the workding dir and submit a new workflow.
        """
        self.clean()
        self.submit()

    def debug(self):
        """
        Initiate an IPython debugging environment to explore the currently
        active SeisFlows environment. Reloads the system modules in an
        interactive environment allowing exploration of the package space.
        Requires 'ipdb' and 'IPython'
        """
        self._load_modules(precheck=False)

        # Distribute modules to common names for easy access during debug mode
        PATH = sys.modules["seisflows_paths"]
        PAR = sys.modules["seisflows_parameters"]
        system = sys.modules["seisflows_system"]
        preprocess = sys.modules["seisflows_preprocess"]
        solver = sys.modules["seisflows_solver"]
        postprocess = sys.modules["seisflows_postprocess"]
        optimize = sys.modules["seisflows_optimize"]
        workflow = sys.modules["seisflows_workflow"]

        # Import debugging options. Following lines will be displayed to stdout
        import ipdb
        from IPython import embed
        # > This is SeisFlows' debug mode.
        ipdb.set_trace(context=5)
        embed(colors="Neutral")
        # > Type 'n' to access a more useful IPython debugger environment.
        # > Type 'workflow.checkpoint()' to save any changes made here.

    def par(self, parameter, value=None):
        """
        Check or set parameters in the SeisFlows parameter file.

        USAGE

            seisflows par [parameter] [value]

            To check the parameter 'NPROC' from the command line:

                seisflows par nproc

            To set the parameter 'BEGIN' to 2:

                seisflows par begin 2

            To change the scratch path to the current working directory:

                seisflows par scratch ./

        :type parameter: str
        :param parameter: parameter to check in parameter file. case insensitive
        :type value: str
        :param value: value to set for parameter. if None, will simply print out
            the current parameter value. To set as Nonetype, set to 'null'
        """
        # SeisFlows parameter file dictates upper-case parameters
        parameter = parameter.upper()

        if not os.path.exists(self._args.parameter_file):
            sys.exit(f"\n\tParameter file '{self._args.parameter_file}' "
                     f"does not exist\n")

        if value is not None and value.lower() == "none":
            warnings.warn("To set values to NoneType, use 'null' not 'None'",
                         UserWarning)

        with open(self._args.parameter_file, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if f"{parameter}: " in line and "#" not in line:
                if value is not None:
                    # These values still have string formatters attached
                    current_par, current_val = line.split(":")

                    # This retains the string formatters of the line
                    new_val = current_val.replace(current_val.strip(), value)
                    lines[i] = ":".join([current_par, new_val])
                    print(f"\n\t{current_par.strip()}: "
                          f"{current_val.strip()} -> {value}\n")

                    with open(self._args.parameter_file, "w") as f:
                        f.writelines(lines)
                else:
                    print(f"\n\t{line}")
                break
        else:
            sys.exit(f"\n\t'{parameter}' not found in parameter file\n")

    def edit(self, name, module, editor=None):
        """
        Directly edit the SeisFlows source code matching the given name
        and module using the chosen text editor.

        USAGE

            seisflows edit [name] [module] [editor]

            To edit the base Solver class using vim, one would run:

                seisflows edit solver base vim

            To simply find the location of the inversion workflow source code:

                seisflows edit workflow inversion q

        :type name: str
        :param name: name of module, must match seisflows.config.NAMES
        :type module: str
        :param module: the module name contained under the SeisFlows namespace
        :type editor: str
        :param editor: optional chosen text editor to open the file.
            * If NoneType: defaults to system environment $EDITOR
            * If 'q': For quit, does not open an editor, simply prints fid
        """
        editor = editor or os.environ.get("EDITOR")
        if editor is None:
            sys.exit("\n\t$EDITOR environment variable not set, set manually\n")

        REPO_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
        if name not in NAMES:
            sys.exit(f"\n\t{name} not in SeisFlows names: {NAMES}\n")

        for package in PACKAGES:
            fid_try = os.path.join(REPO_DIR, package, name, f"{module}.py")
            if os.path.exists(fid_try):
                if editor != "q":
                    subprocess.call([editor, fid_try])
                    sys.exit(f"\n\tEdited file: {fid_try}\n")
                else:
                    sys.exit(f"\n{fid_try}\n")
        else:
            sys.exit(f"\n\tseisflows.{name}.{module} not found\n")

    def check(self, choice=None, *args):
        """
        Check parameters, state or values  of an active SeisFlows environment.
        Type 'seisflows check --help' for a detailed help message.

        :type choice: str
        :param choice: underlying sub-function to choose
        """
        # Initiate a sub parser to provide nested help functions
        subparser = self._parser.add_subparsers(help="subcommand help")
        parser_check = subparser.add_parser(
            name="check", usage=" seisflows check [-h] [subcommand]",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""
Check parameters, state, or values of an active environment
    
    model     check the min/max values of currently active models tracked by
              optimize. 'seisflows check model [name]' to check specific model.
    iter      Check current interation and step count of workflow
    src       List source names and respective internal indices
    isrc      Check source name for corresponding index
            """,
            help="Additional help for the check subcommand")
        parser_check.add_argument("subcommand", type=str, nargs="?",
                                  help="subcommand provided to check command"
                                  )
        parser_check.add_argument("subargs", type=str, default=None, nargs="?",
                                  help="Optional name argument for subcommand")

        if not sys.argv[2:]:
            parser_check.print_help()
            sys.exit(0)
        else:
            # Must ignore the initial positional argument
            args = parser_check.parse_args(sys.argv[2:])

        acceptable_args = {"model": self._check_model_parameters,
                           "iter": self._check_current_iteration,
                           "src": self._check_source_names,
                           "isrc": self._check_source_index}

        if args.subcommand not in acceptable_args.keys():
            parser_check.print_help()
            sys.exit(0)

        self._load_modules(precheck=False)
        acceptable_args[args.subcommand](args.subargs)

    def reset(self, choice=None, *args):
        """
        Mid-level function to wrap lower level reset functions
        """
        acceptable_args = {"line_search": self._reset_line_search,}

        if choice is None:
            sys.exit(f"\n\tseisflows reset requires sub-function from the "
                     f"following:\n\t{list(acceptable_args.keys())}\n")

        if choice not in acceptable_args.keys():
            sys.exit(f"\n\tseisflows reset has no sub-function '{choice}'\n")

        self._load_modules(precheck=False)
        acceptable_args[choice](*args)

    def inspect(self, *args):
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
        self._load_modules(precheck=False)

        if len(args) in [0, 1]:
            self._inspect_module_hierarchy(*args)
        elif len(args) == 2:
            self._inspect_class_that_defined_method(*args)

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
        self._load_modules(precheck=False)

        solver = sys.modules["seisflows_solver"]
        optimize = sys.modules["seisflows_optimize"]
        PATH = sys.modules["seisflows_paths"]

        if path is None:
            path = os.path.join(PATH.OUTPUT, name)
        if os.path.exists(path):
            sys.exit(f"\n\t{path} exists and this action would overwrite the "
                     f"existing path\n")

        solver.save(solver.split(optimize.load(name)), path=path, **kwargs )

    @staticmethod
    def _inspect_class_that_defined_method(name, func):
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
            sys.exit(f"\n\tSeisFlows has no module named '{name}'\n")
        try:
            method = getattr(module, func)
        except AttributeError:
            sys.exit(f"\n\tSeisFlows.{name} has no function '{func}'\n")

        method_name = method.__name__
        if method.__self__:
            classes = [method.__self__.__class__]
        else:
            # Deal with unbound method
            classes = [method.im_class]
        while classes:
            c = classes.pop()
            if method_name in c.__dict__:
                print(f"\n\tSeisFlows.{name}.{func} defined by {c.__name__}\n")
                return
            else:
                classes = list(c.__bases__) + classes
        sys.exit(f"\n\tError matching class for SeisFlows.{name}.{func}\n")

    @staticmethod
    def _inspect_module_hierarchy(name=None):
        """
        Determine the order of class hierarchy for a given SeisFlows module.

        https://stackoverflow.com/questions/1401661/
                            list-all-base-classes-in-a-hierarchy-of-given-class

        :type name: str
        :param name: choice of module, if None, will print hierarchies for all
            modules.
        """
        for name_ in NAMES:
            if name and name_ != name:
                continue
            module = sys.modules[f"seisflows_{name_}"]
            print(f"\n\t{name_.upper()}", end=" ")
            for i, cls in enumerate(inspect.getmro(type(module))[::-1]):
                print(f"-> {cls.__name__}", end=" ")
        print("\n\n\tseisflows inspect [name] [method] to inspect "
              "method ownership\n")

    def _reset_line_search(self):
        """
        Reset the machinery of the line search
        """
        optimize = sys.modules["seisflows_optimize"]
        workflow = sys.modules["seisflows_workflow"]
        
        current_step = optimize.line_search.step_count
        optimize.line_search.reset()
        new_step = optimize.line_search.step_count
    
        print(f"Step Count: {current_step} -> {new_step}")
        workflow.checkpoint()

    def _check_model_parameters(self, src=None):
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
                sys.exit(f"\n\t{src} not in available models {avail}\n")
            srcs = [src]
        for tag in srcs:
            m = optimize.load(tag)
            optimize.check_model_parameters(m, tag)

    def _check_current_iteration(self, *args):
        """
        Display the current point in the workflow in terms of the iteration
        and step count number. Args are not used by allow for a more general
        check() function.
        """
        optimize = sys.modules["seisflows_optimize"]
        try:
            line = optimize.line_search
            cstr = (f"\n"
                    f"\tIteration:  {optimize.iter}\n"
                    f"\tStep Count: {line.step_count} / {line.step_count_max}\n"
                    )
            print(cstr)
        except AttributeError:
            sys.exit("\n\toptimization module has not been initialized yet\n")

    def _check_source_names(self, source_name=None):
        """
        Sources are tagged by name but also by index in the source names which
        can be confusing and usually requires doubling checking. This check
        just prints out source names next to their respective index, or if a
        source name is requested, provides the index for that

        :type source_name: str
        :param source_name: name of source to check index, if None will simply
            print out all sources
        """     
        solver = sys.modules["seisflows_solver"]

        if source_name:
            print(f"{solver.source_names.index(source_name)}: {source_name}")
        else:
            for i, source_name in enumerate(solver.source_names):
                print(f"{i:>3}: {source_name}")

    def _check_source_index(self, idx=None):
        """
        Look up source name by index

        :type idx: int
        :param idx: index of source to look up
        """     
        solver = sys.modules["seisflows_solver"]
        print(f"\n\t{idx}: {solver.source_names[int(idx)]}\n")


def main():
    """
    Main entry point into the SeisFlows package is via the SeisFlows class
    """
    SeisFlows()

