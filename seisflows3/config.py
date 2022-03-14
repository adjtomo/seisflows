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

# List of module names required by SeisFlows for module imports. Order-sensitive
NAMES = ["system", "preprocess", "solver",
         "postprocess", "optimize", "workflow"]

# Packages that define the source code, used to search for base- and subclasses
PACKAGES = ["seisflows3", "seisflows3-super"]

# The location of this config file, which is the main repository
ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# Define a package-wide default directory and file naming schema. This will
# be returned as a Dict() object, defined below.
CFGPATHS = dict(
    PAR_FILE= "parameters.yaml",
    SCRATCHDIR = "scratch",
    STATSDIR = "stats",
    OUTPUTDIR = "output",
    LOGFILE = "output_sf3.txt"
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
    assert("seisflows_parameters" in sys.modules)
    assert("seisflows_paths" in sys.modules)

    # Check if objects already exist on disk
    if os.path.exists(_output()):
        print(msg.WarningOverwrite)
        sys.exit()

    # Instantiate and register objects
    for name in NAMES:
        sys.modules[f"seisflows_{name}"] = custom_import(name)()

    # Error checking
    for name in NAMES:
        sys.modules[f"seisflows_{name}"].check()

    # Ensure that certain parameters are instantiated
    if not hasattr(sys.modules["seisflows_parameters"], "WORKFLOW"):
        print(msg.MissingParameter_Workflow)
        sys.exit(-1)
    if not hasattr(sys.modules["seisflows_parameters"], "SYSTEM"):
        print(msg.MissingParameter_System)
        sys.exit(-1)


def save():
    """
    Export the current session to disk
    """
    logger.info("exporting current working environment to disk")
    unix.mkdir(_output())

    # Save the paths and parameters into a JSON file
    for name in ["parameters", "paths"]:
        fullfile = os.path.join(_output(), f"seisflows_{name}.json")
        savejson(fullfile, sys.modules[f"seisflows_{name}"].__dict__)

    # Save the current workflow as pickle objects
    for name in NAMES:
        fullfile = os.path.join(_output(), f"seisflows_{name}.p")
        saveobj(fullfile, sys.modules[f"seisflows_{name}"])


def load(path):
    """
    Imports a previously saved session from disk

    :type path: str
    :param path: path to the previously saved session
    """
    logger.info("loading current working environment from disk")

    # Load parameters and paths from a JSON file
    for name in ['parameters', 'paths']:
        fullfile = os.path.join(_full(path), f"seisflows_{name}.json")
        sys.modules[f"seisflows_{name}"] = Dict(loadjson(fullfile))

    # Load the saved workflow from pickle objects
    for name in NAMES:
        fullfile = os.path.join(_full(path), f"seisflows_{name}.p")
        sys.modules[f"seisflows_{name}"] = loadobj(fullfile)


class Dict(object):
    """
    Re-defined dictionary-like object for holding parameters or paths

    Allows for easier access of dictionary items, does not allow resets of
    attributes once defined, only allows updates through new dictionaries.
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

    def __init__(self, newdict):
        self.update(newdict)
    
    def __str__(self):
        """
        Pretty print dictionaries and first level nested dictionaries
        """
        str_ = "{"
        for key, item in vars(self).items():
            if isinstance(item, str):
                str_ += f"{key}: '{item}'\n"
            elif isinstance(item, dict):
                str_ += key + ": {\n"
                for key2, item2 in vars(self)[key].items():
                    if isinstance(item2, str):
                        str_ += f"\t{key2}: '{item2}'\n"
                    else:
                        str_ += f"\t{key2}: {item2}\n"
                str_ += "}\n"
            else:
                str_ += f"{key}: {item}\n"
        str_ += "}"
        return str_
        

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
        :param paremeter: name of the parameter
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
            PATH = sys.modules["seisflows_paths"]
            for key, attrs in self.paths.items():
                if attrs["required"] and (key not in PATH):
                    raise ParameterError(PATH, key)
                elif key not in PATH:
                    setattr(PATH, key, attrs["default"])

        if parameters:
            PAR = sys.modules["seisflows_parameters"]
            for key, attrs in self.parameters.items():
                if attrs["required"] and (key not in PAR):
                    raise ParameterError(PAR, key)
                elif key not in PAR:
                    setattr(PAR, key, attrs["default"])


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
        sys.exit(msg.ImportError1)
    # Invalid `system` call
    elif name not in NAMES:
        sys.exit(msg.ImportError2)
    # Attempt to retrieve currently assigned classname from parameters
    if module is None:
        module = _try(name)
        # If no module by that name, return Null
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
    for package in PACKAGES:
        full_dotted_name = ".".join([package, name, module])
        if module_exists(full_dotted_name):
            _exists = True
            break
    if not _exists:
        sys.exit(msg.ImportError3.format(name=name, module=module,
                                         module_upper=name.upper()))
    # Import module
    module = import_module(full_dotted_name)

    # Extract classname from module if possible
    try:
        return getattr(module, classname)
    except AttributeError:
        sys.exit(msg.ImportError4.format(name=name, module=module,
                                         classname=classname))


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


def _par(key):
    """
    Utility function to return a global parameter

    :type key: str
    :param key: key to be passed to the parameter dict
    :return: value for given key
    """
    return sys.modules['seisflows_parameters'][key.upper()]


def _path(key):
    """
    Utility function to return a global pathname

    :type key: str
    :param key: key to be passed to the parameter dict
    :return: value for given key
    """
    return sys.modules['seisflows_paths'][key.upper()]


def _try(key):
    """
    Utility function to try to return a global parameter

    :type key: str
    :param key: key to be passed to the parameter dict
    :return: value for given key or None
    """
    try:
        return _par(key)
    except KeyError:
        return None


def _output():
    """
    Utility function to return the full path to the output directory

    :rtype: str
    :return: full path to output directory
    """
    try:
        return _full(_path("output"))
    except IOError:
        return _full(os.path.join(".", "output"))


def _full(path):
    """
    Utility function to return a full path with trailing '/'

    :type path: str
    :param path: path to be returned in full
    :rtype: str
    :return: full path
    """
    try:
        return os.path.join(os.path.abspath(path), '')
    except Exception as e:
        raise IOError


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
CFGPATHS = Dict(CFGPATHS)
