#!/usr/bin/env python
"""
SeisFlows consists of interacting objects:
'system', 'preprocess', 'solver', 'postprocess', 'optimize', 'workflow'

Each corresponds simultaneously to a module in the SeisFlows source code,
a class that is instantiated and made accessible via sys.modules, and a
parameter in a global dictionary. Once in memory, these objects can be thought
of as comprising the complete 'state' of a SeisFlows session
"""

import os
# import re
import sys
# import imp
import types
import copyreg
from importlib import import_module

from seisflows.tools import msg
# from seisflows.tools.err import ParameterError
from seisflows.tools import unix
from seisflows.tools.tools import loadjson, loadobj, savejson, saveobj
from seisflows.tools.tools import module_exists, package_exists


# The following list is one of the few hardwired aspects of the whole
# SeisFlows package. Any changes may result in circular imports, other problems
names = []
names += ['system']
names += ['preprocess']
names += ['solver']
names += ['postprocess']
names += ['optimize']
names += ['workflow']


def config():
    """
    Instantiates SeisFlows objects and makes them globally accessible by
    registering them in sys.modules
    """
    # Parameters and paths must already be loaded (normally done by sfsubmit)
    assert 'seisflows_parameters' in sys.modules
    assert 'seisflows_paths' in sys.modules

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
        print(msg.MissingParameter_Worfklow)
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
        savejson(fullfile, sys.modules['seisflows_' + name].__dict__)

    # Save the current workflow as pickle objects
    for name in names:
        fullfile = os.path.join(_output(), f"seisflows_{name}.p")
        saveobj(fullfile, sys.modules['seisflows_' + name])


def load():
    """
    Imports a previously saved session from disk
    """
    # Load parameters and paths from a JSON file
    for name in ['parameters', 'paths']:
        fullfile = os.path.join(_output(), f"seisflows_{name}.json")
        sys.modules['seisflows_' + name] = Dict(loadjson(fullfile))

    # Load the saved workflow from pickle objects
    for name in names:
        fullfile = os.path.join(_output(), f"seisflows_{name}.p")
        sys.modules['seisflows_' + name] = loadobj(fullfile)


class Dict(object):
    """
    Re defined dictionary-like object for holding parameters or paths
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


def custom_import(*args):
    """
    Imports SeisFlows module and extracts class of same name.

    For example:
        custom_import('workflow', 'inversion')

        imports 'seisflows.workflow.inversion' and, from this module, extracts
        class 'inversion'.
    """
    # Parse input arguments
    if len(args) == 0:
        raise Exception(msg.ImportError1)
    if args[0] not in names:
        raise Exception(msg.ImportError2)
    if len(args) == 1:
        args += (_try(args[0]),)
    if not args[1]:
        return Null

    # Generate package list
    packages = ['seisflows']

    # Check if modules exist
    _exists = False
    for package in packages:
        full_dotted_name = f"{package}.{args[0]}.{args[1]}"
        if module_exists(full_dotted_name):
            _exists = True
            break
    if not _exists:
        raise Exception(msg.ImportError3 % (args[0], args[1], args[0].upper()))

    # Import module
    module = import_module(full_dotted_name)

    # Extract class
    if hasattr(module, args[1]):
        return getattr(module, args[1])
    else:
        raise Exception(msg.ImportError4 % (args[0], args[1], args[1]))


def tilde_expand(mydict):
    """
    Expands tilde character (~) in Path strings

    :type mydict: dict
    :param mydict: dictionary of paths to be expanded
    """
    for key, val in mydict.items():
        if not isinstance(val, str):
            raise Exception("Expanded objects must type: str")
        elif val[:2] == '~/':
            mydict[key] = f"{os.getenv('HOME')}/{val[2:]}"

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

