#!/usr/bin/env python
"""
The main entry point to the SeisFlows package. A high-level command line tool 
that allows user interface with the underlying SeisFlows environment.
"""
import os
import sys
import argparse
from glob import glob

from seisflows.tools import unix, tools
from seisflows.config import config, tilde_expand, Dict, names
from seisflows.tools.tools import loadyaml


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
    parser.add_argument("-p", "--parameter_file", nargs="?",
                        default="parameters.yaml")

    return parser.parse_args()


def parse_null(parameters):
    """
    Seisflows was written in such a way that excluding parameters from the 
    parameter file would enforce default values. This becomes confusing however 
    because then the User needs to look at the definitions in the check()
    functions of each module to determine which parameters they want to exclude.

    Rather than rewrite this system, we allow the User to set parameters to 
    None, null, '', etc. This function sanitizes those inputs before handing 
    them over to Seisflows. This allows for full parameter files that include 
    all options, but that also works within the existing machinery of SeisFlows.

    This function also enforces a few string formats 

    :type parameters: dict
    :param parameters: dict of parameters to parse
    :rtype: dict
    :return: parameters that have been sanitized of all null values
    """
    # Copy the dictionary to get around deleting keys while iterating
    parsed_parameters = dict(parameters)
    for key, item in parameters.items():
        # Search for all None and "" items, ignore bools, 0's etc.
        if not item and isinstance(item, (type(None), str)):
            del parsed_parameters[key]
        elif key == "LINESEARCH":
            parsed_parameters[key] = item.capitalize()

    return parsed_parameters


def setup(precheck=True):
    """
    Load the paths and parameters from file into sys modules, set the default
    parameters if they are missing from the file, and expand all paths to 
    absolute pathnames.

    :type precheck: bool
    :param precheck: print out a few key parameters and require user-input 
        before allowing workflow to be submitted. This is usually run before
        submit and resume, to prevent job submission without evaluation.
    """
    args = get_args()

    # Check if the filepaths exist
    if not os.path.exists(args.parameter_file):
        raise FileNotFoundError(f"Parameter file not found: "
                                f"{args.parameter_file}")

    # Register parameters
    parameters = loadyaml(args.parameter_file)
    if precheck:
        precheck_parameters(parameters)
    parameters = parse_null(parameters)
    sys.modules["seisflows_parameters"] = Dict(parameters)

    # Register paths, expand to relative paths to absolute
    paths = tilde_expand(parameters["PATHS"])
    paths = {key: os.path.abspath(path) for key, path in paths.items()}
    sys.modules["seisflows_paths"] = Dict(paths)

    return args, paths, parameters


def load_modules(**kwargs):
    """
    A function to load and check each of the SeisFlows modules, initiating
    the SeisFlows environment. All modules are reliant on one another so any
    access to SeisFlows requires loading everything simultaenously.
    Calls setup() beforehand which loads in paths and parameters as well.
    """
    args, paths, parameters = setup(**kwargs)
    
    # Working directory should already have been created by submit()
    unix.cd(args.workdir)

    # Reload objects from Pickle files
    for name in names:
        fullfile = os.path.join(args.workdir, "output",
                                f"seisflows_{name}.p")
        sys.modules[f"seisflows_{name}"] = tools.loadobj(fullfile)

    # Check parameters so that default values are present
    for name in names:
        sys.modules[f"seisflows_{name}"].check()


def precheck_parameters(par):
    """
    Print important arguments before beginning a workflow. This will ensure that
    the User is aware of how key arguments are set that will have major affect
    on the workflow.
    """
    msg = f"""
      Beginning workflow "{par['TITLE']}"...
          WORKFLOW: {par['WORKFLOW']}
          BEGIN: {par['BEGIN']}
          END: {par['END']}
          RESUME_FROM: {par['RESUME_FROM']}
          STOP_AFTER: {par['STOP_AFTER']}
          NTASK: {par['NTASK']}
          WALLTIME: {par['WALLTIME']}
          CASE: {par['CASE']}
          SMOOTH H,V: {par['SMOOTH_H']}, {par['SMOOTH_V']}
          """
    print(msg)
    check = input("Continue? (y/[n]): ")
    if check != "y":
        sys.exit(-1)


def submit():
    """
    Submit the workflow for the first time. Create the working directory and 
    any required paths and ensure that all required paths exist.
    """
    args, paths, parameters = setup(precheck=True)

    # Check that paths exist
    paths_dont_exist = []
    for key, path in paths.items():
        if (key in ["OUTPUT", "SCRATCH"]) or (path == ""):
            continue
        if not os.path.exists(path):
            paths_dont_exist.append(path)
    if paths_dont_exist:
        print("\nThe following paths do not exist:\n")
        for path_ in paths_dont_exist:
            print(f"\t{path_}")
        print("\n")
        sys.exit()

    unix.mkdir(args.workdir)
    unix.cd(args.workdir)

    # Submit workflow
    config()  # Instantiate all modules and check their parameters
    workflow = sys.modules["seisflows_workflow"]
    system = sys.modules["seisflows_system"]
    system.submit(workflow)


def clean():
    """
    Clean the working directory by deleting everything except the parameter file
    """
    args = get_args() 
    check = input("\nThis will remove all workflow objects, leaving only the "
                  "parameter file.\nAre you sure you want to clean? (y/[n]): ")
    if check == "y":
        for fid in glob(os.path.join(args.workdir, "output*")):
            unix.rm(fid)
        for fid in glob(os.path.join(args.workdir, "*log*")):
            unix.rm(fid)
        unix.rm(os.path.join(args.workdir, "scratch"))


def resume():
    """
    Resume a previously started workflow by loading the module pickle files and
    submitting the workflow from where it left off.
    """
    load_modules(precheck=True)

    workflow = sys.modules["seisflows_workflow"]
    system = sys.modules["seisflows_system"]

    system.submit(workflow)


def restart():
    """
    Restart simply means clean the cwd and submit a new workflow
    """
    clean()
    submit()


def debug():
    """
    A debug mode that reloads the system modules and starts an IPython
    debugger allowing exploration of the package space in an interactive 
    environment.
    """
    load_modules(precheck=False)
    
    # Distribute modules to common names for easy access during debug mode
    PATH = sys.modules["seisflows_paths"]
    PAR = sys.modules["seisflows_parameters"]
    system = sys.modules["seisflows_system"]
    preprocess = sys.modules["seisflows_preprocess"]
    solver = sys.modules["seisflows_solver"]
    postprocess = sys.modules["seisflows_postprocess"]
    optimize = sys.modules["seisflows_optimize"]
    workflow = sys.modules["seisflows_workflow"]

    # Import debugging options. The following lines will be displayed to console
    import ipdb
    from IPython import embed
    # This is Seisflows' debug mode.
    ipdb.set_trace(context=5)
    embed(colors="Neutral")
    # type 'n' to access a more useful IPython debugger.
    # type 'workflow.checkpoint()' to save any changes made here.


class SeisShows:
    """
    A high-level API that allows the user to manipulate or inspect the SeisFlows 
    workflow via simple command line arguments. Almost every modules requires 
    loading of other modules so to run checks we must load the entire SeisFlows
    environment, so this is slower than it could be, but it provides the most 
    flexibility when accessing internal information.

    Throws sys.exit if invalid arguments are passed to the 'seisflows' command
    """
    def __init__(self):
        """
        Initiate the entire SeisFlows workflow and module system then call the
        requested function provided through command line arguments
        """ 
        args = get_args()

        # If SeisShows is being called, then the format for calling it splits
        # the argument into (func, *args). Dont care about the other argparse
        try:
            func, *extra_args = args.main_args

            if not hasattr(self, func):
                # Available functions are any public methods
                sys.exit(f"\n\tseisflows has no matching function '{func}'"
                         f"\n\tbesides main functions, available sub-functions "
                         f"include:\n\t{self._public_methods}\n")

            getattr(self, func)(*extra_args)
        # If class initiated with no arguments, assume its for debug reasons
        except ValueError:
            pass

    @property
    def _public_methods(self):
        """Return a list of all public methods within this class"""
        return [_ for _ in dir(self) if not _.startswith("_")]

    def check(self, choice=None, *args):
        """
        Mid-level function to call more lower level check functions 
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

        load_modules()
        acceptable_args[choice](*args)

    def reset(self, choice=None, *args):
        """
        Mid-level function to call lower level reset functions
        """
        acceptable_args = {"line_search": self._reset_line_search,}

        if choice is None:
            sys.exit(f"\n\tseisflows reset requires sub-function from the "
                     f"following:\n\t{list(acceptable_args.keys())}\n")

        if choice not in acceptable_args.keys():
            sys.exit(f"\n\tseisflows reset has no sub-function '{choice}'\n")

        load_modules()
        acceptable_args[choice](*args)
        
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
    Main entry point into the SeisFlows package
    """
    # Easy way to convert strings to functions
    acceptable_args = {"submit": submit, "resume": resume, "clean": clean,
                       "restart": restart, "debug": debug}

    try:
        main_arg = get_args().main_args[0]
    except IndexError:
        ok_args = list(acceptable_args.keys()) + SeisShows()._public_methods
        sys.exit(f"\n\tseisflows command requires an argument. "
                 f"available arguments are:\n\t{ok_args}\n"
                 f"\ttype 'seisflows -h' for a help message\n")

    if main_arg in acceptable_args.keys():
        acceptable_args[main_arg]()
    else:
        SeisShows()


if __name__ == "__main__":
    main()
        
