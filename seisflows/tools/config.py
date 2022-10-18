#!/usr/bin/env python3
"""
Seisflows configuration tools, containing core utilities that are called upon
throughout the Seisflows workflow.
"""
import dill
import logging
import os
import sys
import re
import yaml
import numpy as np
import traceback
from pkgutil import find_loader
from importlib import import_module

from seisflows import logger, NAMES
from seisflows.tools import msg

ENV_VARIABLES = ["SEISFLOWS_TASKID", "SLURM_ARRAY_TASK_ID"]


class Dict(dict):
    """
    A dictionary replacement which allows for easier parameter access through
    getting and setting attributes. Also has some functionality to make string
    printing prettier
    """
    def __str__(self):
        """Pretty print dictionaries and first level nested dictionaries"""
        str_ = ""
        try:
            longest_key = max([len(_) for _ in self.keys()])
            for key, val in self.items():
                str_ += f"{key:<{longest_key}}: {val}\n"
        except ValueError:
            pass
        return str_

    def __repr__(self):
        """Pretty print when calling an instance of this object"""
        return self.__str__()

    def __getattr__(self, key):
        """Attribute-like access of the internal dictionary attributes"""
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"{key} not found in Dict")

    def __setattr__(self, key, val):
        """Setting attributes can only be performed one time"""
        self.__dict__[key] = val


class Null:
    """
    A null object that always and reliably does nothing
    """
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        return self

    def __bool__(self):
        return False

    def __nonzero__(self):
        return False

    def __getattr__(self, key):
        return self

    def __setattr__(self, key, val):
        return self

    def __delattr__(self, key):
        return self


def load_yaml(filename):
    """
    Define how the PyYaml yaml loading function behaves.
    Replaces None and inf strings with NoneType and numpy.inf respectively

    :type filename: str
    :param filename: .yaml file to load in
    :rtype: Dict
    :return: Dictionary containing all parameters in a YAML file
    """
    # work around PyYAML bugs
    yaml.SafeLoader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.'))

    with open(filename, 'r') as f:
        mydict = Dict(yaml.safe_load(f))

    if mydict is None:
        mydict = Dict()

    # Replace 'None' and 'inf' values to match expectations
    for key, val in mydict.items():
        if val == "None":
            mydict[key] = None
        if val == "inf":
            mydict[key] = np.inf

    return mydict


def get_task_id():
    """
    Task IDs are assigned to each child process spawned by the system module
    during a SeisFlows workflow. SeisFlows modules use this Task ID to keep
    track of embarassingly parallel process, e.g., solver uses the Task ID to
    determine which source is being considered.

    :rtype: int
    :return: task id for given solver
    """
    for env_var in ENV_VARIABLES:
        _taskid = os.getenv(env_var)
        if _taskid is not None:
            return int(_taskid)
    else:
        logger.warning("Environment Task ID variable not found. Assigning 0")
        return 0


def set_task_id(task_id):
    """
    Set the SEISFLOWS_TASKID in os environs for local workflows. If running
    on HPC systems, running array jobs will assign the Task ID

    .. note::
        Mostly used for debugging/testing purposes as a way of mimicing
        system.run() assigning task ids to child processes

    :type task_id: int
    :param task_id: integer task id to assign to the current working environment
    """
    os.environ["SEISFLOWS_TASKID"] = str(task_id)


def import_seisflows(workdir=os.getcwd(), parameter_file="parameters.yaml",
                     **kwargs):
    """
    Standard SeisFlows workflow setup block which runs a number of setup
    tasks including: loading a user-defined parameter file, configuring the
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
                  verbose=parameters.verbose, **kwargs)

    # Instantiate SeisFlows modules dynamically based on choices and parameters
    # provided in the input parameter file
    modules = Dict()
    for name in NAMES[:]:
        # Workflow is instantiated differently
        if name == "workflow":
            continue
        modules[name] = custom_import(name, parameters[name])(**parameters)

    # Import workflow separately by providing all the instantiated modules to it
    workflow = \
        custom_import("workflow", parameters["workflow"])(modules, **parameters)

    return workflow


def config_logger(level="DEBUG", filename=None, filemode="a", verbose=True,
                  stream_handler=True):
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

    # Stream handler to print log statements to stdout. Sometimes we don't want
    # this, e.g., on an HPC system having both stream and file will print
    # double log messages to your log file
    if stream_handler:
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


def pickle_function_list(functions, path=os.getcwd(), **kwargs):
    """
    Save a list of functions and their keyword arguments as pickle files.
    Return the names of the files. Used for running functions from spawned
    processes during cluster runs.

    .. note::
        The idea here is that we need this list of functions to be
        discoverable by a system separate to the one that defined them. To
        do this we can pickle Python objects on disk, and have the new
        system read in the pickle files and evaluate the objects. We use
        'dill' because Pickle can't serialize methods/functions

    :type functions: list of methods
    :param functions: a list of functions that should be run in order. All
        kwargs passed to run() will be passed into the functions.
    :type path: str
    :param path: path to save the pickle files. Defaults to current working
        directory
    :rtype: tuple of str
    :return: (name of the pickle file containing the function,
        name of the pickle file containing keyword arguments)
    """
    # Save the instances that define the functions as a pickle object
    func_names = "_".join([_.__name__ for _ in functions])  # unique identifier
    fid_funcs_pickle = os.path.join(path, f"{func_names}.p")

    with open(fid_funcs_pickle, "wb") as f:
        dill.dump(obj=functions, file=f)

    # Save the kwargs as a separate pickle object
    fid_kwargs_pickle = os.path.join(path, f"{func_names}_kwargs.p")
    with open(fid_kwargs_pickle, "wb") as f:
        dill.dump(obj=kwargs, file=f)

    return fid_funcs_pickle, fid_kwargs_pickle


def number_fid(fid, i=0):
    """
    Number a filename. Used to store old log files without overwriting them.
    Premise is, if you have a file e.g., called: output.txt
    This function would return incrementing filenames:
    output_000.txt, output_001.txt, output_002.txt, ouput_003.txt ...

    .. note::
        Replace statement is catch-all, so we assume that there is only one
        instance of the file extension in the entire path.

    :type fid: str
    :param fid: path to the file that you want to increment
    :type i: int
    :param i: number to append to file id
    :rtype: str
    :return: filename with appended number. filename ONLY, will strip away
        the original path location
    """
    fid_only = os.path.basename(fid)
    ext = os.path.splitext(fid_only)[-1]  # e.g., .txt
    new_ext = f"_{i:0>3}{ext}"   # e.g., _000.txt
    new_fid = fid_only.replace(ext, new_ext)
    return new_fid
