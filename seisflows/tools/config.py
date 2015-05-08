import os
import imp
import sys
import types
from os.path import abspath, join

from seisflows.tools import unix
from seisflows.tools.code import Struct, loadobj, savejson, saveobj


class ConfigObj(object):
    """ SeisFlows consists of interacting 'system', 'preprocess', 'solver',
      'postprocess', 'optimize', and 'workflow' objects. The role of the
      ConfigObj utility is to initialize these objects and make them accessible
      via the standard Python import system.

      Objects are first created via the 'ConfigObj.initialize' method.  In the 
      course of executing a workflow, it is often necessary to read or write 
      objects to disk, for example, to copy an object to a remote node or to 
      pause and resume executation. Methods 'ConfigObj.save' and 
      'ConfigObj.load' are provided for this purpose.

      Objects themselves are customizable, with different choices available for
      each object. For example, in the main package, two workflows, 'inversion'
      and 'migration', are provided. Users can choose between the two through 
      the 'WORKFLOW' setting in the  parameters file. If desired functionality
      is missing from the main package, users can create their own workflows or
      customize existing workflows.
"""
    # Because each 'object' is itself customizable, it is rarely necessary to 
    # modify the following list. If it modifications are desired anyway, caution
    # should be excercised, as changing the names of objects or the order in
    # which they loaded can result in circular imports or other problems.
    objects = []
    objects += ['system']
    objects += ['preprocess']
    objects += ['solver']
    objects += ['postprocess']
    objects += ['optimize']
    objects += ['workflow']

    def __init__(self, name='SeisflowsObjects'):
        self.name = name

    def __iter__(self):
        return iter(list(self.objects))

    def init(self):
        """ Instantiates objects, registers objects in sys.modules, and calls 
            objects' check methods
        """
        try:
            parameters = sys.modules['SeisflowsParameters']
        except:
            raise Exception

        # instantiate objects
        for obj in self.objects:
            par = parameters[obj.upper()]
            sys.modules[obj] = loadclass(obj, par)()

        # check objects
        for obj in self.objects:
            sys.modules[obj].check()

    def register(self, key, val):
        """ Makes object globally accessible by registering it in sys.modules
        """
        #if key in sys.modules:
        #    raise Exception

        sys.modules[key] = val

    def unregister(self, key):
        """ Removes object from sys.modules
        """
        sys.modules.pop(key)

    def save(self, name, path='.'):
        """ Saves all objects to disk
        """
        try:
            fullpath = join(abspath(path), name)
        except:
            raise IOError(path)

        unix.mkdir(fullpath)
        for key in self.objects:
            saveobj(fullpath + '/' + key+'.p', sys.modules[key])

    def load(self, name, path='.'):
        """ Loads objects from disk
        """
        try:
            fullpath = join(abspath(path), name)
        except:
            raise IOError(path)

        for obj in self.objects:
            fullfile = join(fullpath, obj+'.p')
            sys.modules[obj] = loadobj(fullfile)

        for obj in self.objects:
            sys.modules[obj].check()



class ParameterObj(object):
    """ Dictionary like object for holding parameters. Makes parameters globally 
        accessible by registering itself in sys.modules
    """

    def __new__(self, name, path=None):
        if name in sys.modules:
            return sys.modules[name]
        else:
            return object.__new__(self)

    def __init__(self, name, path=None):
        if name not in sys.modules:
            sys.modules[name] = self

    def init(self, path):
        self.update(loadvars(path, '.'))
        return self

    def __iter__(self):
        return iter(sorted(self.__dict__.keys()))

    def __getattr__(self, key):
        return self.__dict__[key]

    #def __getitem__(self, key):
    #    return self.__dict__[key]

    def __setattr__(self, key, val):
        if key in self.__dict__:
            raise TypeError("Once defined, parameters cannot be changed.")
        self.__dict__[key] = val

    def __delattr__(self, key):
        if key in self.__dict__:
            raise TypeError("Once defined, parameters cannot deleted.")
        raise KeyError

    def update(self, newdict):
        super(ParameterObj, self).__setattr__('__dict__', newdict)

    def save(self, name):
        savejson(name, self.__dict__)


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


class ParameterError(ValueError):
    def __init__(self, obj, key):
        if key not in obj:
            message = '%s is not defined.' % key
        elif key in obj:
            message = '%s has bad value: ' % key, obj.__getattr__(key)

        super(ParameterError, self).__init__(message)


def loadclass(*args):
    """ Given name of module relative to package directory, returns
        corresponding class

        Note: The module should have a class with the exact same name.
    """
    if not args:
        return Null

    if not args[-1]:
        return Null

    # first, try importing relative to main package directory
    list = _parse(args, package='seisflows')
    string = '.'.join(list)
    if _exists(list):
        obj = getattr(_import(string), list[-1])
        return obj

    # next, try importing relative to extensions directory
    list = _parse(args, package='seisflows_research')
    string = '.'.join(list)
    if _exists(list):
        obj = getattr(_import(string), list[-1])
        return obj

    raise ImportError('Not found in SeisFlows path.')


def loadvars(*args, **kwargs):
    return _vars(_import(*args, **kwargs))


def findpath(obj):
    """ Determines absolute path of
            - a module, from an instance
            - a file, from its name
            - a seisflow module, from its full (dotted) name
    """
    if isinstance(obj, types.ModuleType):
        path = obj.__file__
    elif os.path.isfile(obj):
        path = os.path.abspath(obj)
    else:
        string = '.'.join(_parse([obj], package='seisflows'))
        moduleobj = _import(string)
        path = moduleobj.__file__
    return os.path.dirname(path)


# -- utility functions

def _import(string, path=None):
    """Imports from string"""
    if path:
        # temporarily adjust python path
        sys.path.append(path)

    moduleobj = __import__(string, fromlist='dummy')

    if path:
        sys.path.pop()

    return moduleobj


def _vars(obj):
    """Returns an object's __dict__ with private variables removed"""
    mydict = {}
    for key, val in vars(obj).items():
        if key[0] != '_':
            mydict[key] = val
    return Struct(mydict)


def _parse(args, package=None):
    """Determines path of module relative to package directory"""
    arglist = []
    for arg in args:
        arglist.extend(arg.split('.'))
    if package:
        parglist = package.split('.')
        nn = min(len(arglist), len(parglist))
        for ii in range(nn + 1):
            if parglist[ii:nn] == arglist[0:nn - ii]:
                break
        arglist = parglist[0:ii] + arglist
    return arglist


def _exists(parts):
    """Checks if a module exists without importing it

    Takes an array of strings.
    """
    try:
        path = None
        for part in parts[:-1]:
            args = imp.find_module(part, path)
            obj = imp.load_module(part, *args)
            path = obj.__path__
        args = imp.find_module(parts[-1], path)
        return True
    except ImportError:
        return False

