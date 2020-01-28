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

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.tools import loadjson, loadobj, savejson, saveobj
from seisflows.tools.tools import module_exists


# The following list is one of the few hardwired aspects of the whole
# SeisFlows package. Any changes may result in circular imports, other problems
names = ["system", "preprocess", "solver", "postprocess",
         "optimize", "workflow"]


def config():
    """
    Instantiates SeisFlows objects and makes them globally accessible by
    registering them in sys.modules
    """
    # Parameters and paths must already be loaded (normally done by submit)
    assert('seisflows_parameters' in sys.modules)
    assert('seisflows_paths' in sys.modules)

    # Check if objects already exist on disk
    if os.path.exists(_output()):
        print(msg.WarningOverwrite)
        sys.exit()

    # Instantiate and register objects
    for name in names:
        sys.modules['seisflows_' + name] = custom_import(name)()

    # Error checking
    for name in names:
        sys.modules['seisflows_' + name].check()

    # Ensure that certain parameters are instantiated
    if not hasattr(sys.modules['seisflows_parameters'], "WORKFLOW"):
        print(msg.MissingParameter_Workflow)
        sys.exit(-1)
    if not hasattr(sys.modules['seisflows_parameters'], "SYSTEM"):
        print(msg.MissingParameter_System)
        sys.exit(-1)


def save():
    """
    Export the current session to disk
    """
    unix.mkdir(_output())

    # Save the paths and parameters into a JSON file
    for name in ['parameters', 'paths']:
        fullfile = os.path.join(_output(), f"seisflows_{name}.json")
        savejson(fullfile, sys.modules[f"seisflows_{name}"].__dict__)

    # Save the current workflow as pickle objects
    for name in names:
        fullfile = os.path.join(_output(), f"seisflows_{name}.p")
        saveobj(fullfile, sys.modules[f"seisflows_{name}"])


def load(path):
    """
    Imports a previously saved session from disk

    :type path: str
    :param path: path to the previously saved session
    """
    # Load parameters and paths from a JSON file
    for name in ['parameters', 'paths']:
        fullfile = os.path.join(_full(path), f"seisflows_{name}.json")
        sys.modules[f"seisflows_{name}"] = Dict(loadjson(fullfile))

    # Load the saved workflow from pickle objects
    for name in names:
        fullfile = os.path.join(_full(path), f"seisflows_{name}.p")
        sys.modules[f"seisflows_{name}"] = loadobj(fullfile)


class Dict(object):
    """
    Re defined dictionary-like object for holding parameters or paths

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
                str_ += (f"{key}: '{item}'\n")
            elif isinstance(item, dict):
                str_ += (key + ": {\n")
                for key2, item2 in vars(self)[key].items():
                    if isinstance(item2, str):
                        str_ += (f"\t{key2}: '{item2}'\n")
                    else:
                        str_ += (f"\t{key2}: {item2}\n")
                str_ += "}\n"
            else:
                str_ += (f"{key}: {item}\n")
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


def custom_import(name=None, module=None, classname=None):
    """
    Imports SeisFlows module and extracts class that is the camelcase version
    of the module name

    For example:
        custom_import('workflow', 'inversion')

        imports 'seisflows.workflow.inversion' and, from this module, extracts
        class 'inversion'.

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
        raise Exception(msg.ImportError1)
    # Invalid `system` call
    elif name not in names:
        raise Exception(msg.ImportError2)
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

    # Generate package list
    packages = ["seisflows"]

    # Check if modules exist, otherwise raise custom exception
    _exists = False
    for package in packages:
        full_dotted_name = ".".join([package, name, module])
        if module_exists(full_dotted_name):
            _exists = True
            break
    if not _exists:
        raise Exception(
            msg.ImportError3.format(name=name, module=module,
                                    module_upper=module.upper())
        )
    # Import module
    module = import_module(full_dotted_name)

    # Extract classname from module if possible
    try:
        return getattr(module, classname)
    except AttributeError:
        raise Exception(
            msg.ImportError4.format(name=name, module=module,
                                    classname=classname)
        )


def tilde_expand(mydict):
    """
    Expands tilde character (~) in Path strings

    :type mydict: dict
    :param mydict: dictionary of paths to be expanded
    """
    for key, val in mydict.items():
        if not isinstance(val, str):
            raise Exception("Expanded objects must type: str")
        else:
            mydict[key] = os.path.expanduser(val)

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

