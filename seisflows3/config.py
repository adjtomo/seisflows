#!/usr/bin/env python
"""
This is the Seisflows Config script, it contains utilities that are called upon
throughout the Seisflows workflow. It also (re)defines some important functions
that are used extensively by the machinery of Seisflows.

SeisFlows consists of interacting objects:
'system', 'preprocess', 'solver', 'postprocess', 'optimize', 'workflow'

Each corresponds simultaneously to a module in the SeisFlows source code,
a class that is instantiated and made accessible via sys.modules, and a
parameter in a global dictionary. Once in memory, these objects can be thought
of as comprising the complete 'state' of a SeisFlows session
"""
import os
import sys
import types
import copyreg
from importlib import import_module

from seisflows3 import logger
from seisflows3.tools import msg, unix
from seisflows3.tools.wrappers import loadjson, loadobj, savejson, saveobj
from seisflows3.tools.wrappers import module_exists
from seisflows3.tools.err import ParameterError


"""
!!! WARNING !!!

The following constants are (some of the only) hardwired components
of the pacakge. The naming, order, case, etc., of each constant may be 
important, and any changes to these will more-than-likely break the underlying 
mechanics of the package. Do not touch unless you know what you're doing!
"""

# List of module names required by SeisFlows3 for imports. Order-sensitive
# In sys.modules these will be prepended by 'seisflows_', e.g., seisflows_system
NAMES = ["system", "preprocess", "solver",
         "postprocess", "optimize", "workflow"]

# Packages that define the source code, used to search for base- and subclasses
PACKAGES = ["seisflows3", "seisflows3-super"]

# These define the sys.modules names where parameter values and paths are stored
PAR = "seisflows_parameters"
PATH = "seisflows_paths"

# The location of this config file, which is the main repository
ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# Define a package-wide default directory and file naming schema. This will
# be returned as a Dict() object, defined below. All of these files and
# directories will be created relative to the user-defined working directory
CFGPATHS = dict(
    PAR_FILE="parameters.yaml",  # Default SeisFlows3 parameter file
    SCRATCHDIR="scratch",        # SeisFlows3 internal working directory
    STATSDIR="stats",            # Optimization module log file output
    OUTPUTDIR="output",          # Permanent disk storage for state and outputs
    LOGFILE="output_sf3.txt",    # Log files for all system log
    ERRLOGFILE="error_sf3.txt",  # StdErr dump site for crash messages
    LOGDIR="logs",               # Dump site for previously created log files
)
"""
!!! ^^^ WARNING ^^^ !!!
"""


def init_seisflows():
    """
    Instantiates SeisFlows3 objects and makes them globally accessible by
    registering them in sys.modules
    """
    logger.info("initializing SeisFlows3 in sys.modules")

    # Parameters and paths must already be loaded (normally done by submit)
    assert(PAR in sys.modules)
    assert(PATH in sys.modules)

    # Check if objects already exist on disk, exit so as to not overwrite
    if "OUTPUT" in sys.modules[PATH] and \
            os.path.exists(sys.modules[PATH]["OUTPUT"]):
        print(msg.cli("Data from previous workflow found in working directory.",
                      items=["> seisflows restart: delete data and start new "
                             "workflow",
                             "> seisflows resume: resume existing workflow"],
                      header="warning", border="=")
              )
        sys.exit(-1)

    # Instantiate and register objects
    for name in NAMES:
        sys.modules[f"seisflows_{name}"] = custom_import(name)()

    # Error checking
    for name in NAMES:
        sys.modules[f"seisflows_{name}"].check()

    # Bare minimum module requirements for SeisFlows3
    req_modules = ["WORKFLOW", "SYSTEM"]
    for req in req_modules:
        if not hasattr(sys.modules[PAR], req):
            print(msg.cli(f"SeisFlows3 requires defining: {req_modules}."
                          "Please specify these in the parameter file. Use "
                          "'seisflows print module' to determine suitable "
                          "choices.", header="error", border="="))
            sys.exit(-1)


def save():
    """
    Export the current session to disk
    """
    logger.info("exporting current working environment to disk")
    output = sys.modules[PATH]["OUTPUT"]
    unix.mkdir(output)

    # Save the paths and parameters into a JSON file
    for name in [PAR, PATH]:
        fullfile = os.path.join(output, f"{name}.json")
        savejson(fullfile, sys.modules[name].__dict__)

    # Save the current workflow as pickle objects
    for name in NAMES:
        fullfile = os.path.join(output, f"seisflows_{name}.p")
        saveobj(fullfile, sys.modules[f"seisflows_{name}"])


def load(path):
    """
    Imports a previously saved session from disk

    :type path: str
    :param path: path to the previously saved session
    """
    logger.info("loading current working environment from disk")

    # Load parameters and paths from a JSON file
    for name in [PAR, PATH]:
        fullfile = os.path.join(os.path.abspath(path), f"{name}.json")
        sys.modules[name] = Dict(loadjson(fullfile))

    # Load the saved workflow from pickle objects
    for name in NAMES:
        fullfile = os.path.join(os.path.abspath(path), f"seisflows_{name}.p")
        sys.modules[f"seisflows_{name}"] = loadobj(fullfile)


class Dict(object):
    """
    Re-defined dictionary-like object for holding parameters or paths

    Allows for easier access of dictionary items, does not allow resets of
    attributes once defined, only allows updates through new dictionaries.

    !!! TODO Not sure if it makes sense to define an entirely new object here
    !!! TODO as we are essentially re-inventing the wheel? Why not inhereit
    !!! TODO dict, i.e., Dict(dict), and then modify a few of the key attributes
    !!! TODO that we don't want, like setattr and delattr. Dict is pretty
    !!! TODO central to seisflows though so I don't want to break anything
    !!! TODO without sufficient tests in place to figure out what broke.
    """
    def __iter__(self):
        return iter(sorted(self.__dict__.keys()))

    def __getattr__(self, key):
        return self.__dict__[key]

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setattr__(self, key, val):
        if key in self.__dict__:
            raise TypeError("Once defined, parameters cannot be changed.")
        self.__dict__[key] = val

    def __delattr__(self, key):
        if key in self.__dict__:
            raise TypeError("Once defined, parameters cannot be deleted.")
        raise KeyError

    def update(self, newdict):
        super(Dict, self).__setattr__('__dict__', newdict)

    def values(self):
        return self.__dict__.values()

    def __init__(self, newdict):
        self.update(newdict)
    
    def __str__(self):
        """
        Pretty print dictionaries and first level nested dictionaries
        """
        str_ = ""
        longest_key = max([len(_) for _ in self.__dict__.keys()])
        for key, val in self.__dict__.items():
            str_ += f"{key:<{longest_key}}: {val}\n"
        return str_

    def __repr__(self):
        """
        Pretty print
        :return:
        """
        return(self.__str__())
        

class Null(object):
    """
    A null object that always and reliably does nothing
    """
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        return self

    def __nonzero__(self):
        return False

    def __getattr__(self, key):
        return self

    def __setattr__(self, key, val):
        return self

    def __delattr__(self, key):
        return self


class SeisFlowsPathsParameters:
    """
    A class used to simplify defining required or optional paths and parameters
    by enforcing a specific structure to their entry into the environment.
    Replaces the functionalities of the old check() functions.

    .. note::
        if a path or parameter is optional it requires a default value.
    """
    def __init__(self, base=None):
        """
        We simply store paths and parameters as nested dictioanries. Due to the
        use of inheritance, the class can be passed to itself on initialization
        which means paths and parameters can be adopted from base class

        :type base: seisflows.config.DefinePathsParameters
        :param base: paths and parameters from abstract Base class that need to 
            be inherited by the current child class.
        """
        self.parameters, self.paths = {}, {}
        if base:
            self.parameters.update(base.parameters)
            self.paths.update(base.paths)

    def par(self, parameter, required, docstr, par_type, default=None):
        """
        Add a parameter to the internal list of parameters

        :type parameter: str
        :param parameter: name of the parameter
        :type required: bool
        :param required: whether or not the parameter is required. If it is not
            required, then a default value should be given
        :type docstr: str
        :param docstr: Short explanatory doc string that defines what the
            parameter is used for.
        :type par_type: class or str
        :param par_type: the parameter type, used for doc strings and also
            parameter validation
        :param default: default value for the parameter, can be any type
        """
        if required:
            default = "!!! REQUIRED PARAMETER !!!"
        if type(par_type) == type:
            par_type = par_type.__name__
        self.parameters[parameter] = {"docstr": docstr, "required": required,
                                      "default": default, "type": par_type}

    def path(self, path, required, docstr, default=None):
        """
        Add a path to the internal list of paths

        :type path: str
        :param path: name of the parameter
        :type required: bool
        :param required: whether or not the path is required. If it is not
            required, then a default value should be given
        :type docstr: str
        :param docstr: Short explanatory doc string that defines what the
            path is used for.
        :type default: str
        :param default: default value for the path

        """
        if required:
            default = "!!! REQUIRED PATH !!!"
        self.paths[path] = {"docstr": docstr, "required": required,
                            "default": default}

    def validate(self, paths=True, parameters=True):
        """
        Set internal paths and parameter values into sys.modules. Should be
        called by each modules check() function.

        Ensures that required paths and parameters are set by the user, and that
        default values are stored for any optional paths and parameters.

        :type paths: bool
        :param paths: validate the internal path values
        :type parameters: bool
        :param parameters: validate the internal parameter values
        :raises ParameterError: if a required path or parameter is not set by
            the user.
        """
        if paths:
            sys_path = sys.modules[PATH]
            for key, attrs in self.paths.items():
                if attrs["required"] and (key not in sys_path):
                    raise ParameterError(sys_path, key)
                elif key not in sys_path:
                    setattr(sys_path, key, attrs["default"])

        if parameters:
            sys_par = sys.modules[PAR]
            for key, attrs in self.parameters.items():
                if attrs["required"] and (key not in sys_par):
                    raise ParameterError(sys_par, key)
                elif key not in sys_par:
                    setattr(sys_par, key, attrs["default"])


def custom_import(name=None, module=None, classname=None):
    """
    Imports SeisFlows module and extracts class that is the camelcase version
    of the module name

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
            "name 'seisflows3.name.module' must correspond to a module "
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
            module = sys.modules[PAR][name.upper()]
        except KeyError:
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
    for package in PACKAGES:
        full_dotted_name = ".".join([package, name, module])
        if module_exists(full_dotted_name):
            _exists = True
            break
    if not _exists:
        print(msg.cli(f"The following module was not found within the package: "
                      f"seisflows3.{name}.{module}",
                      header="custom import error", border="=")
              )
        sys.exit(-1)
    # Import module
    module = import_module(full_dotted_name)

    # Extract classname from module if possible
    try:
        return getattr(module, classname)
    except AttributeError:
        print(msg.cli(f"The following method was not found in the imported "
                      f"class: seisflows3.{name}.{module}.{classname}"))
        sys.exit(-1)


def format_paths(mydict):
    """
    Ensure that paths have a standardized format before being allowed into
    an active working environment. 
    Expands tilde character (~) in path strings and expands absolute paths

    :type mydict: dict
    :param mydict: dictionary of paths to be expanded
    :rtype: dict
    :return: formatted path dictionary
    """
    for key, val in mydict.items():
        try:
            mydict[key] = os.path.expanduser(os.path.abspath(val))
        except TypeError:
            continue
    return mydict


def _pickle_method(method):
    """
    The following code changes how instance methods are handled by pickle.
    Placing it in this module ensures that pickle changes will be in
    effect for all SeisFlows workflows

    Note: For relevant discussion, see stackoverflow thread:
    "Can't pickle <type 'instancemethod'> when using python's
    multiprocessing Pool.map()"

    Relevant Links (last accessed 01.20.2020):
        https://stackoverflow.com/questions/7016567/
        picklingerror-when-using-multiprocessing

        https://bytes.com/topic/python/answers/
        552476-why-cant-you-pickle-instancemethods
    """
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    """
    The unpickling counterpart to the above function
    """
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)


copyreg.pickle(types.MethodType, _pickle_method, _unpickle_method)

# Because we defined Dict inside this file, we need to convert our CFGPATHS
# to a Dict at the end of the file to allow direct variable access
# !!! TODO I don't really like this implementation, can it be changed?
CFGPATHS = Dict(CFGPATHS)
