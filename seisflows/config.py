#!/usr/bin/env python3
"""
This is the Seisflows Config script, it contains utilities that are called upon
throughout the Seisflows workflow. It also (re)defines some important functions
that are used extensively by the machinery of Seisflows.


Each corresponds simultaneously to a module in the SeisFlows source code,
a class that is instantiated and made accessible via sys.modules, and a
parameter in a global dictionary. Once in memory, these objects can be thought
of as comprising the complete 'state' of a SeisFlows session
"""
import os
import sys
import logging
import traceback
from pkgutil import find_loader
from importlib import import_module

from seisflows import logger
from seisflows.tools.core import Dict, Null
from seisflows.tools import msg
from seisflows.tools.core import load_yaml


# List of module names required by SeisFlows for imports. Order-sensitive
NAMES = ["workflow", "system", "solver", "preprocess", "optimize"]

# The location of this config file, which is the main repository
ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def import_seisflows(workdir=os.getcwd(), parameter_file="parameters.yaml"):
    """
    Standard SeisFlows workflow setup block which runs a number of setup
    tasks including: loading a user-defiend parameter file, configuring the
    package-wide logger based on user-input path to log file and desired
    verbosity, and instantiating all modules in a generic fashion based on user
    choice. Returns the 'workflow' module, which contains all other submodules
    as attributes.

    :type workdir: str
    :param workdir: the current working directory in which to perform a
        SeisFlows workflow. Defaults to the current working directory
    :type parameter_file: str
    :param parameter_file: the YAML formatted parameter file that is used to
        instantiate each of the SeisFlows modules and run the workflow. This
        should be created by the command line argument 'seisflows configure'.
        Defaults to 'parameters.yaml'
    :rtype: module
    :return: instantiated 'workflow' module which contains all sub-modules which
        have been instantiated with user-defined parameters
    """
    # Read in parameters from file. Set up the logger
    parameters = load_yaml(os.path.join(workdir, parameter_file))
    config_logger(level=parameters.log_level,
                  filename=parameters.path_output_log,
                  verbose=parameters.verbose)

    # Instantiate SeisFlows modules dynamically based on choices and parameters
    # provided in the input parameter file
    modules = Dict()
    for name in NAMES[:]:
        # Workflow is instantiated differently
        if name == "workflow":
            continue
        modules[name] = custom_import(name, parameters[name])(**parameters)
        # parameters.pop(name)  # drop name so workflow doesnt instantiate it

    # Import workflow separately by providing all the instantiated modules to it
    workflow = \
        custom_import("workflow", parameters["workflow"])(modules, **parameters)

    return workflow


def config_logger(level="DEBUG", filename=None, filemode="a", verbose=True):
    """
    Explicitely configure the logging module with some parameters defined
    by the user in the System module. Instantiates a stream logger to write
    to stdout, and a file logger which writes to `filename`. Two levels of
    verbosity and three levels of log messages allow the user to determine
    how much output they want to see.

    :type level: str
    :param level: log level to be passed to logger, available are
        'CRITICAL', 'WARNING', 'INFO', 'DEBUG'
    :type filename: str or None
    :param filename: name of the log file to write log statements to. If None,
        logs will be written to STDOUT ONLY, and `filemode` will not be used.
    :type filemode: str
    :param filemode: method for opening the log file. defaults to append 'a'
    :type verbose: bool
    :param verbose: if True, writes a more detailed log message stating the
        type of log (warning, info, debug), and the class and method which
        called the logger (e.g., seisflows.solver.specfem2d.save()). This
        is much more useful for debugging but clutters up the log file.
        if False, only write the time and message in the log statement.
    """
    # Make sure that we don't already have handlers described, which may happen
    # if this function gets run multiple times, and leads to duplicate logs
    while logger.hasHandlers() and logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Two levels of verbosity on log level, triggered with PAR.VERBOSE
    if verbose:
        # More verbose logging statement for debugging
        fmt_str = (
            "%(asctime)s | %(levelname)-5s | "
            "%(filename)s -> %(funcName)s():L%(lineno)s\n"
            "> %(message)s"
        )
    else:
        # Clean logging statement with only time and message
        fmt_str = "%(asctime)s (%(levelname).1s) | %(message)s"

    # Instantiate logger during _register() as we now have user-defined pars
    logger.setLevel(level)
    formatter = logging.Formatter(fmt_str, datefmt="%Y-%m-%d %H:%M:%S")

    # Stream handler to print log statements to stdout
    st_handler = logging.StreamHandler(sys.stdout)
    st_handler.setFormatter(formatter)
    logger.addHandler(st_handler)

    # File handler to print log statements to text file `filename`
    if filename is not None:
        file_handler = logging.FileHandler(filename, filemode)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)


def custom_import(name=None, module=None, classname=None):
    """
    Imports SeisFlows module and extracts class that is the camelcase version
    of the module name. Used to dynamically import sub-modules by name only,
    avoiding the need to hardcode import statements.

    For example:
        custom_import('workflow', 'inversion')

        imports 'seisflows.workflow.inversion' and, from this module, extracts
        class 'Inversion'.

    :type name: str
    :param name: component of the workflow to import, defined by `names`,
        available: "system", "preprocess", "solver",
                   "postprocess", "optimize", "workflow"
    :type module: module within the workflow component to call upon, e.g.
        seisflows.workflow.inversion, where `inversion` is the module
    :type classname: str
    :param classname: the class to be called from the module. Usually this is
        just the CamelCase version of the module, which will be defaulted to if
        this parameter is set `None`, however allows for custom class naming.
        Note: CamelCase class names following PEP-8 convention.
    """
    # Parse input arguments for custom import
    # Allow empty system to be called so that import error message can be thrown
    if name is None:
        print(msg.cli(
            "Please check that 'custom_import' utility is being used as "
            "follows: custom_import(name, module). The resulting full dotted "
            "name 'seisflows.name.module' must correspond to a module "
            "within this package.", header="custom import error", border="="))
        sys.exit(-1)
    # Invalid `system` call
    elif name not in NAMES:
        print(msg.cli(
            "Please check that the use of custom_import(name, module, class) "
            "is implemented correctly, where name must be in the following:",
            items=NAMES, header="custom import error", border="="))
        sys.exit(-1)
    # Attempt to retrieve currently assigned classname from parameters
    if module is None:
        try:
            module = sys.modules["seisflows_parameters"][name.upper()]
        except KeyError:
            return Null
        # If this still returns nothing, then no module has been assigned
        # likely the User has turned this module OFF
        if module is None:
            return Null
    # If no method specified, convert classname to PEP-8
    if classname is None:
        # Make a distinction for fully uppercase classnames, e.g. LBFGS
        if module.isupper():
            classname = module.upper()
        # If normal classname, convert to CamelCase
        else:
            classname = module.title().replace("_", "")

    # Check if modules exist, otherwise raise custom exception
    _exists = False
    full_dotted_name = ".".join(["seisflows", name, module])
    # find_loader() checks if the module exists or not
    if not find_loader(full_dotted_name):
        print(msg.cli(f"The following module was not found within the package: "
                      f"seisflows.{name}.{module}",
                      header="custom import error", border="=")
              )
        sys.exit(-1)

    # If importing the module doesn't work, throw an error. Usually this happens
    # when an external dependency isn't available, e.g., Pyatoa
    try:
        module = import_module(full_dotted_name)
    except Exception as e:
        print(msg.cli(f"Module could not be imported {full_dotted_name}",
                      items=[str(e)], header="custom import error", border="="))
        print(traceback.print_exc())
        sys.exit(-1)

    # Extract classname from module if possible
    try:
        return getattr(module, classname)
    except AttributeError:
        print(msg.cli(f"The following method was not found in the imported "
                      f"class: seisflows.{name}.{module}.{classname}"))
        sys.exit(-1)

# def save(path):
#     """
#     Export the current Python environment to disk as Pickle and JSON files,
#     which allows us to checkpoint a current workflow and resume without
#     loss of information.
#
#     :type path: str
#     :param path: path to save the current session
#     """
#     if not os.path.exists(path):
#         unix.mkdir(path)
#
#     # Save the paths and parameters into a JSON file
#     for name in ["seisflows_parameters", "seisflows_paths"]:
#         fullfile = os.path.join(path, f"{name}.json")
#         with open(fullfile, "w") as f:
#             json.dump(sys.modules[name], f, sort_keys=True, indent=4)
#
#     # Save the current workflow as pickle objects
#     for name in NAMES:
#         fullfile = os.path.join(path, f"seisflows_{name}.p")
#         with open(fullfile, "wb") as f:
#             pickle.dump(sys.modules[f"seisflows_{name}"], f)
#
#
# def load(path):
#     """
#     Imports a previously saved session from disk by reading in JSON and
#     Pickle files which define a saved Python environment
#
#     :type path: str
#     :param path: path to the previously saved session
#     """
#     # Load parameters and paths from a JSON file
#     for name in ["seisflows_parameters", "seisflows_paths"]:
#         fullfile = os.path.join(os.path.abspath(path), f"{name}.json")
#         with open(fullfile, "r") as f:
#             sys.modules[name] = Dict(json.load(f))
#
#     # Load the saved workflow from pickle objects
#     for name in NAMES:
#         fullfile = os.path.join(os.path.abspath(path), f"seisflows_{name}.p")
#         with open(fullfile, "rb") as f:
#             sys.modules[f"seisflows_{name}"] = pickle.load(f)

# def _pickle_method(method):
#     """
#     The following code changes how instance methods are handled by pickle.
#     Placing it in this module ensures that pickle changes will be in
#     effect for all SeisFlows workflows
#
#     Note: For relevant discussion, see stackoverflow thread:
#     "Can't pickle <type 'instancemethod'> when using python's
#     multiprocessing Pool.map()"
#
#     Relevant Links (last accessed 01.20.2020):
#         https://stackoverflow.com/questions/7016567/
#         picklingerror-when-using-multiprocessing
#
#         https://bytes.com/topic/python/answers/
#         552476-why-cant-you-pickle-instancemethods
#     """
#     func_name = method.im_func.__name__
#     obj = method.im_self
#     cls = method.im_class
#     return _unpickle_method, (func_name, obj, cls)
#
#
# def _unpickle_method(func_name, obj, cls):
#     """
#     The unpickling counterpart to the above function
#     """
#     for cls in cls.mro():
#         try:
#             func = cls.__dict__[func_name]
#         except KeyError:
#             pass
#         else:
#             break
#     return func.__get__(obj, cls)


# copyreg.pickle(types.MethodType, _pickle_method, _unpickle_method)

