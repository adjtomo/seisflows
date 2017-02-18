
import copy_reg
import imp
import os
import re
import sys
import types

from importlib import import_module
from pkgutil import find_loader
from os.path import abspath, join, exists

from seisflows.tools import msg
from seisflows.tools.err import ParameterError
from seisflows.tools import unix
from seisflows.tools.tools import Struct, loadjson, loadobj, loadpy, savejson, saveobj

# SeisFlows consists of interacting 'system', 'preprocess', 'solver', 'postprocess', 'optimize', and 'workflow' objects. Each corresponds simultaneously to a module in the SeisFlows source code, a class that is instantiated and made accessible via sys.modules, and a parameter in a global dictionary. Once in memory, these objects can be thought of as comprising the complete 'state' of a SeisFlows session

# The following list is one of the few hardwired aspects of the whole SeisFlows package. Any changes may result in circular imports or other problems

names = []
names += ['system']
names += ['preprocess']
names += ['solver']
names += ['postprocess']
names += ['optimize']
names += ['workflow']


def config():
    """ Instantiates SeisFlows objects and makes them globally accessible by
      registering them in sys.modules
   """
    # check if objects exist on disk
    if exists(_full(_path())):
        print msg.WarningOverwrite
        sys.exit()

    if 'seisflows_parameters' not in sys.modules:
        raise Exception

    if 'seisflows_paths' not in sys.modules:
        raise Exception

    # instantiate and register objects
    for name in names:
        sys.modules['seisflows_'+name] = custom_import(name)()

    # parameter checking
    for name in names:
        sys.modules['seisflows_'+name].check()


def save():
    """ Exports session to disk
    """
    unix.mkdir(_full(_path()))

    for name in ['parameters', 'paths']:
        fullfile = join(_full(_path()), 'seisflows_'+name+'.json')
        savejson(fullfile, sys.modules['seisflows_'+name].__dict__)

    for name in names:
        fullfile = join(_full(_path()), 'seisflows_'+name+'.p')
        saveobj(fullfile, sys.modules['seisflows_'+name])


def load(path):
    """ Imports session from disk
    """
    for name in ['parameters', 'paths']:
        fullfile = join(_full(path), 'seisflows_'+name+'.json')
        sys.modules['seisflows_'+name] = Dict(loadjson(fullfile))

    for name in names:
        fullfile = join(_full(path), 'seisflows_'+name+'.p')
        sys.modules['seisflows_'+name] = loadobj(fullfile)


class Dict(object):
    """ Dictionary-like object for holding parameters or paths
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
    """ Always and reliably does nothing
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
    """ Imports SeisFlows module and extracts class of same name. For example,

            custom_import('workflow', 'inversion') 

        imports 'seisflows.workflow.inversion' and, from this module, extracts
        class 'inversion'.
    """
    # parse input arguments
    if len(args) == 0:
        raise Exception(msg.ImportError1)
    if args[0] not in names:
        raise Exception(msg.ImportError2)
    if len(args) == 1:
        args += (_val(args[0]),)
    if not args[1]:
        return Null

    # generate package list
    packages = ['seisflows']
    if os.getenv('SEISFLOWS_PACKAGES'):
        for package in os.getenv('SEISFLOWS_PACKAGES').split(','):
            if package in packages:
                continue
            if find_loader(package):
                packages += [package]

    # does module exist?
    _exists = False
    for package in packages:
        full_dotted_name = package+'.'+args[0]+'.'+args[1]
        if find_loader(full_dotted_name):
            _exists = True
            break
    if not _exists:
        raise Exception(msg.ImportError3 % 
            (args[0], args[1], args[0].upper()))

    # import module
    module = import_module(full_dotted_name)

    # extract class
    if hasattr(module, args[1]):
        return getattr(module, args[1])
    else:
        raise Exception(msg.ImportError4 % 
            (args[0], args[1], args[1]))


def tilde_expand(mydict):
    """ Expands tilde character in path strings
    """
    for key,val in mydict.items():
        if type(val) not in [str, unicode]:
            raise Exception
        if val[0:2] == '~/':
            mydict[key] = os.getenv('HOME') +'/'+ val[2:]
    return mydict



# utility functions

def _path():
    try:
        return sys.modules['seisflows_paths']['OUTPUT']
    except:
        cwd = abspath('.')
        return join(cwd, 'output')


def _full(path):
    try:
        return join(abspath(path), '')
    except:
        raise IOError


def _val(key):
    try:
        return sys.modules['seisflows_parameters'][key.upper()]
    except KeyError:
        return None


# the following code changes how instance methods are handled by pickle.  placing it here, in this module, ensures that pickle changes will be in effect for all SeisFlows workflows

# for relevant discussion, see stackoverflow thread "Can't pickle <type 'instancemethod'> when using python's multiprocessing Pool.map()"

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)


copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

