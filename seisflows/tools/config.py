
import copy_reg
import imp
import os
import sys
import types

from os.path import abspath, join

from seisflows.tools import unix
from seisflows.tools.code import Struct, loadjson, loadobj, savejson, saveobj
from seisflows.tools.msg import WarningOverwrite


class SeisflowsObjects(object):
    """ SeisFlows consists of interacting 'system', 'preprocess', 'solver',
      'postprocess', 'optimize', and 'workflow' objects. The role of the
      SeisflowsObjects utility is to load these objects into memory and make 
      them accessible via the standard Python import system. Once in memory,
      these objects can be thought of as comprising the complete 'state' of a
      SeisFlows session.

      Each entry in 'objects' below corresponds simultaneously to a parameter in
      SeisflowsParameters; a module in the SeisFlows package; and a class that
      is instantiated and made accessible via the Python import system.

      The list 'objects' is is one of the few hardwired aspects of the whole
      SeisFlows package. Any changes may result in circular imports or other
      problems.
    """

    objects = []
    objects += ['system']
    objects += ['preprocess']
    objects += ['solver']
    objects += ['postprocess']
    objects += ['optimize']
    objects += ['workflow']


    def load(self):
        """ Load objects from the above list into memeory
        """
        if 'SeisflowsParameters' not in sys.modules.keys():
            raise Exception

        if 'SeisflowsPaths' not in sys.modules.keys():
            raise Exception

        # check if objects from previous run exist on disk
        if os.path.isdir(self.fullpath()):
            print WarningOverwrite
            sys.exit()

        for obj in self.objects:
            # instantiate objects
            key = sys.modules['SeisflowsParameters'][obj.upper()]

            # make them accessible via the stanard Python import system
            sys.modules[obj] = loadclass(obj, key)()

        self.check()


    def save(self, path):
        """ Save objects to disk for later reference
        """
        fullpath = self.fullpath(path)
        unix.mkdir(fullpath)
        for key in self.objects:
            saveobj(fullpath +'/'+ key+'.p', sys.modules[key])


    def reload(self, path):
        """ Load saved objects from disk
        """
        fullpath = self.fullpath(path)
        for obj in self.objects:
            fullfile = join(fullpath, obj+'.p')
            sys.modules[obj] = loadobj(fullfile)
        self.check()


    def check(self):
        for obj in ['system', 'optimize', 'workflow', 'solver', 'preprocess', 'postprocess']:
        #for obj in self.objects:
            sys.modules[obj].check()


    def fullpath(self, path=None):
        if not path:
            try:
                path = sys.modules['SeisflowsParameters']['OUTPUT']
            except:
                path = './output'

        try:
            fullpath = join(abspath(path), 'SeisflowsObjects')
        except:
            raise IOError(path)

        return fullpath


class ParameterObj(object):
    """ Dictionary like object for holding parameters
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
            raise TypeError("Once defined, parameters cannot deleted.")
        raise KeyError

    def update(self, newdict):
        super(ParameterObj, self).__setattr__('__dict__', newdict)


class SeisflowsParameters(ParameterObj):
    def __new__(self):
        if 'SeisflowsParameters' in sys.modules:
            return sys.modules['SeisflowsParameters']
        else:
            return object.__new__(self)

    def __init__(self):
        if 'SeisflowsParameters' not in sys.modules:
            sys.modules['SeisflowsParameters'] = self

    def load(self):
        mydict = loadvars('parameters', '.')
        self.update(mydict)
        return self

    def save(self, path):
        fullfile = join(path, 'SeisflowsParameters.json')
        savejson(fullfile, self.__dict__)

    def reload(self, path):
        fullfile = join(path, 'SeisflowsParameters.json')
        mydict = loadjson(fullfile)
        self.update(mydict)


class SeisflowsPaths(ParameterObj):
    def __new__(self):
        if 'SeisflowsPaths' in sys.modules:
            return sys.modules['SeisflowsPaths']
        else:
            return object.__new__(self)

    def __init__(self):
        if 'SeisflowsPaths' not in sys.modules:
            sys.modules['SeisflowsPaths'] = self

    def load(self):
        mydict = loadvars('paths', '.')
        super(ParameterObj, self).__setattr__('__dict__', mydict)
        return self

    def save(self, path):
        fullfile = join(path, 'SeisflowsPaths.json')
        savejson(fullfile, self.__dict__)

    def reload(self, path):
        fullfile = join(path, 'SeisflowsPaths.json')
        mydict = loadjson(fullfile)
        self.update(mydict)


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
    def __init__(self, *args):
        if len(args) == 0:
            msg = 'Bad parameter.'
            super(ParameterError, self).__init__(msg)
        elif len(args) == 1:
            msg = 'Bad parameter: %s' % args[0]
            super(ParameterError, self).__init__(msg)
        elif args[0] not in args[1]:
            msg = '%s is not defined.' % args[0]
            super(ParameterError, self).__init__(msg)
        elif key in obj:
            msg = '%s has bad value: ' % args[0], args[1].__getattr__(args[0])
            super(ParameterError, self).__init__(msg)


def loadclass(*args):
    """ Given name of module relative to package directory, returns
        corresponding class

        Note: The module should have a class with exactly the same name.
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


### utility functions

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


# the following changes how instance methods are handled by pickle.  placing it here, in this module, ensures that pickle changes will be in effect for all SeisFlows workflows

copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

