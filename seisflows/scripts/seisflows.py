#!/usr/bin/env python
"""
The main entry point to the SeisFlows package. A high-level command line tool 
that allows user interface with the underlying SeisFlows environment.
"""
import os
import sys
import inspect
import argparse
from glob import glob
from textwrap import wrap
from seisflows.tools import unix, tools
from seisflows.tools.tools import loadyaml, loadpy
from seisflows.config import (init_seisflows, tilde_expand, Dict, names, 
                              custom_import, ROOT_DIR, PAR_FILE)


def get_args():
    """
    Get User defined arguments, or assign defaults

    :rtype: argparse.ArgumentParser()
    :return: User defined or default arguments
    """
    parser = argparse.ArgumentParser()

    # Positional arguments
    parser.add_argument("main_args", type=str, nargs="*",
                        help="task for Seisflows to perform")

    # Optional parameters
    parser.add_argument("-w", "--workdir", nargs="?", default=os.getcwd())
    parser.add_argument("--path_file", nargs="?", default="paths.py")
    parser.add_argument("-p", "--parameter_file", nargs="?",
                        default="parameters.yaml")

    return parser.parse_args()


def _parse_null(dictionary):
    """
    Remove null, None or '' values from a dictionary

    :type dictionary: dict
    :param dictionary: dict of parameters to parse
    :rtype: dict
    :return: dictionary that has been sanitized of all null values
    """
    # Copy the dictionary to get around deleting keys while iterating
    parsed_dict = dict(dictionary)
    for key, item in dictionary.items():
        # Search for all None and "" items, ignore bools, 0's etc.
        if not item and isinstance(item, (type(None), str)):
            del parsed_dict[key]

    return parsed_dict


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
        Parse user-defined arguments
        """
        self._args = get_args()

        # To be set by _register()
        self._paths = None
        self._parameters = None

        # If SeisFlows is being called, then the format for calling it splits
        # the argument into (func, *args). Ignore the other argparse arguments
        try:
            func, *extra_args = self._args.main_args

            if not hasattr(self, func):
                # Available functions are any public methods
                sys.exit(f"\n\tseisflows has no matching function '{func}'"
                         f"\n\ttype 'seisflows' for available functions\n")

            getattr(self, func)(*extra_args)
        # If class initiated with no arguments, assume its for debug reasons
        except ValueError:
            pass

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
            raise FileNotFoundError(f"Parameter file not found: "
                                    f"{self._args.parameter_file}")

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
            import warnings
            warnings.warn(".py parameter and path files are deprecated in "
                          "favor of a .yaml parameter file. Please consider "
                          "switching as the use of legacy .py files may have "
                          "unintended consequences at runtime",
                          DeprecationWarning)

            assert(os.path.exists(self._args.paths_file)), \
                f"Legacy parameter file requires corresponding path file"
            parameters = loadpy(self._args.parameter_file)
            paths = loadpy(self._args.path_file)
        else:
            raise TypeError(f"Unknown file format for "
                            f"{self._args.parameter_file}, file must be '.yaml' "
                            f"(preferred) or '.py' (legacy)")

        # WORKDIR needs to be set here as it's expected by most modules
        if "WORKDIR" not in paths:
            paths["WORKDIR"] = self._args.workdir

        # For submit() and resume(), provide a dialogue and require a visual
        # pre-check before submitting the workflow
        if precheck:
            msg = f"""
              Beginning workflow "{parameters['TITLE']}"...
                  WORKFLOW: {parameters['WORKFLOW']}
                  BEGIN: {parameters['BEGIN']}
                  END: {parameters['END']}
                  RESUME_FROM: {parameters['RESUME_FROM']}
                  STOP_AFTER: {parameters['STOP_AFTER']}
                  NTASK: {parameters['NTASK']}
                  WALLTIME: {parameters['WALLTIME']}
                  CASE: {parameters['CASE']}
                  SMOOTH H,V: {parameters['SMOOTH_H']}, {parameters['SMOOTH_V']}
                  """
            print(msg)
            check = input("Continue? (y/[n]): ")
            if check != "y":
                sys.exit(-1)

        # Register parameters to sys, ensure they meet standards of the package
        # parameters = _parse_null(parameters)
        sys.modules["seisflows_parameters"] = Dict(parameters)

        # Register paths to sys, expand to relative paths to absolute, drop null
        paths = tilde_expand(_parse_null(paths))
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
        for name in names:
            fullfile = os.path.join(self._args.workdir, "output",
                                    f"seisflows_{name}.p")
            sys.modules[f"seisflows_{name}"] = tools.loadobj(fullfile)

        # Check parameters so that default values are present
        for name in names:
            sys.modules[f"seisflows_{name}"].check()

    def setup(self):
        """
        Setup a SeisFlows working directory by copying in a template parameter
        file and symlinking the source code for easy access to the repo.
        """
        if os.path.exists(self._args.parameter_file):
            print(f"\n\tParameter file '{self._args.parameter_file}' "
                  f"already exists\n")
            check = input(f"\tOverwrite with blank file? (y/[n]): ")
            if check == "y":
                unix.rm(self._args.parameter_file)
            else:
                sys.exit()

        # Template parameter file should be located in the main directory
        unix.cp(PAR_FILE, self._args.workdir)

        # Symlink the source code for easy access to repo
        if not os.path.exists(os.path.join(self._args.workdir, "source_code")):
            unix.ln(ROOT_DIR, os.path.join(self._args.workdir, "source_code"))

    def configure(self):
        """
        Dynamically generate the parameter file by writing out docstrings and
        default values for each of the SeisFlows modules. Writes the files
        manually, but consistent with the .yaml file format.
        """
        self._register(precheck=False)

        # Set some re-usable strings to provide a consistent look
        BLANK = "\n"
        COMMENT = "#\n"
        DIVIDER = f"# {'=' * 78}\n"
        HEADER = "# {} {}\n"
        UNDERLINE = f"# {'-' * 25}\n"
        TAB = "    "

        HEADER_TOP = BLANK + DIVIDER + COMMENT + HEADER + UNDERLINE + COMMENT
        HEADER_BOT = COMMENT + DIVIDER

        # Establish the paths and parameters provided by the user
        assert(self._args.parameter_file.endswith(".yaml")), \
            f"seisflows configure only applicable to .yaml parameter files"

        # Need to attempt importing all modules before we access any of them
        for name in names:
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
                for name in names:
                    req = sys.modules[f"seisflows_{name}"].required
                    seisflows_paths.update(req.paths)

                    # Write the docstring header with all parameters
                    f.write(HEADER_TOP.format(name.upper(), "PARAMETERS"))
                    for key, attrs in req.parameters.items():
                        f.write(f"# {key} ({attrs['type']}):\n")
                        # Ensure that header lines are no more than 80 char
                        docstrs = wrap(attrs["docstr"], width=80-len(TAB),
                                       break_long_words=False)
                        for line, docstr in enumerate(docstrs):
                            f.write(f"#{TAB}{docstr}\n")
                    f.write(HEADER_BOT)
                    # Write parameters in a YAML consistent format
                    for key, attrs in req.parameters.items():
                        # Yaml writes None as 'null'
                        if attrs["default"] is None:
                            default = "null"
                        else:
                            default = attrs["default"]
                        f.write(f"{key}: {default}\n")

                # Write the paths in the same format as parameters
                f.write(HEADER_TOP.format("PATHS", ""))
                for key, attrs in seisflows_paths.items():
                    # Ensure that header lines are no more than 80 char
                    docstr_ = f"{key}: {attrs['docstr']}"
                    docstrs = wrap(docstr_, width=79 - len(TAB),
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
        Initiate a SeisFlows environment and save the initial state in pickles.
        Similar to submit() but with no path checking and no workflow
        submission. Useful for exploring the environment prior to a workflow
        submission or for debug and development purposes.
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
        Submit the workflow for the first time. Create the working directory and
        any required paths and ensure that all required paths exist.
        """
        self._register(precheck=True)

        # A list of paths that need to exist if provided by user
        REQ_PATHS = ["SPECFEM_BIN", "SPECFEM_DATA", "MODEL_INIT", "MODEL_TRUE",
                     "DATA", "LOCAL", "MASK"]

        # Check that paths exist
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

        # Submit workflow
        init_seisflows()
        workflow = sys.modules["seisflows_workflow"]
        system = sys.modules["seisflows_system"]

        system.submit(workflow)

    def clean(self):
        """
        Clean the working directory by deleting everything except parameter file
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
        Restart simply means clean the workding dir and submit a new workflow
        """
        self.clean()
        self.submit()

    def debug(self):
        """
        A debug mode that reloads the system modules and starts an IPython
        debugger allowing exploration of the package space in an interactive
        debugging environment.
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
        # > Type 'n' to access a more useful IPython debugger.
        # > Type 'workflow.checkpoint()' to save any changes made here.

    def check(self, choice=None, *args):
        """
        Mid-level function to wrap more lower level check functions
        """
        acceptable_args = {"model": self._check_model_parameters,
                           "iter": self._check_current_iteration,
                           "src": self._check_source_names,
                           "isrc": self._check_source_index}

        if choice is None:
            sys.exit(f"\n\tseisflows check requires sub-function from the "
                     f"following:\n\t{list(acceptable_args.keys())}\n")

        if choice not in acceptable_args.keys():
            sys.exit(f"\n\tseisflows check has no sub-function '{choice}'\n")

        self._load_modules(precheck=False)
        acceptable_args[choice](*args)

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
        Inspect inheritance hierarchy of the package for easier debugging
        """
        self._load_modules(precheck=False)

        if len(args) in [0, 1]:
            self._inspect_module_hierarchy(*args)
        elif len(args) == 2:
            self._inspect_class_that_defined_method(*args)

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
        for name_ in names:
            if name and name_ != name:
                continue
            module = sys.modules[f"seisflows_{name_}"]
            print(f"\n{name_.upper()}")
            for i, cls in enumerate(inspect.getmro(type(module))[::-1]):
                print(f"{(i + 1) * '  '}{cls.__name__}")

    # @staticmethod
    # def _inspect_module_hierarchy(name=None):
    #     """
    #     This doesn't work, only prints out the new methods but doesnt tell
    #     if methods have been overwritten
    #     """
    #     for name_ in names:
    #         if name and name_ != name:
    #             continue
    #         module = sys.modules[f"seisflows_{name_}"]
    #         print(f"\n{name_.upper()}")
    #         # Start from the type and move to Baseclass and each Subclass
    #         method_list_sub, method_list_super = [], []
    #         for cls in inspect.getmro(type(module))[::-1]:
    #             # Get a list of public methods from the class that do not match
    #             # the parent class methods, in order to define a list of unique
    #             # methods only defined in this subclass
    #             method_list_super = [func for func in dir(cls) if
    #                                  callable(getattr(cls, func))
    #                                  and not func.startswith("_")
    #                                  and not func in method_list_sub]
    #             print(f"{'=' * 80}\n{cls.__name__}\n{'=' * 80}")
    #             for method in method_list_super:
    #                 print(method)
    #             # Keep track of already defined methods
    #             method_list_sub += method_list_super

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
            assert(src in srcs), f"{src} not in available models {avail}"
            srcs = [src]
        for tag in srcs:
            m = optimize.load(tag)
            optimize.check_model_parameters(m, tag)

    def _check_current_iteration(self):
        """
        Display the current point in the workflow in terms of the iteration
        and step count number
        """
        optimize = sys.modules["seisflows_optimize"]

        line = optimize.line_search
        cstr = (f"\n"
                f"\tIteration:  {optimize.iter}\n"
                f"\tStep Count: {line.step_count} / {line.step_count_max}\n")
        print(cstr)

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
    Main entry point into the SeisFlows package via the SeisFlows class
    """
    try:
        # Force help message if no arguments provided with command
        get_args().main_args[0]
    except IndexError:
        acceptable_args = "\n\t\t".join(SeisFlows()._public_methods)
        sys.exit(f"\n\tseisflows command requires an argument. "
                 f"available arguments are:\n"
                 f"\n\t\t{acceptable_args}\n"
                 f"\n\ttype 'seisflows -h' for a help message\n")

    SeisFlows()


if __name__ == "__main__":
    main()
        
