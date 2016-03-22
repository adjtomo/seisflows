
import copy_reg
import imp
import os
import re
import sys
import types

from importlib import import_module
from pkgutil import find_loader
from os.path import abspath, join, exists

from seisflows.tools import unix
from seisflows.tools.code import Struct, loadjson, loadobj, savejson, saveobj, loadpy
from seisflows.tools.msg import WarningOverwrite, ImportError1, ImportError2, ImportError3, ImportError4


class SeisflowsObjects(object):
    """ SeisFlows consists of interacting 'system', 'preprocess', 'solver',
      'postprocess', 'optimize', and 'workflow' objects. The role of the
      SeisflowsObjects utility is to instantiate these objects and make them
      accessible via the standard Python import system. Once in memory, these
      objects should be thought of as comprising the complete 'state' of a
      SeisFlows session.

      The following list is one of the few hardwired aspects of the whole 
      SeisFlows package. Any changes may result in circular imports or other
      problems. Each entry corresponds simultaneously to a parameter in
      SeisflowsParameters; a module in the SeisFlows package; and a class that
      is instantiated and made accessible via the Python import system.
    """
    names = []
    names += ['system']
    names += ['optimize']
    names += ['preprocess']
    names += ['solver']
    names += ['postprocess']
    names += ['workflow']


    def load(self):
        """ Imports and instantiates objects from the above list
        """
        if 'SeisflowsParameters' not in sys.modules.keys():
            raise Exception

        if 'SeisflowsPaths' not in sys.modules.keys():
            raise Exception

        # check if objects from previous run exist on disk
        if exists(_full(_path())):
            print WarningOverwrite
            sys.exit()

        for name in self.names:
            # import and instantiate object
            obj = custom_import(name)()

            # make accessible via standard Python import system
            sys.modules[name] = obj

        self.check()


    def save(self, path):
        """ Save objects to disk for later reference
        """
        unix.mkdir(_full(path))
        for name in self.names:
            fullfile = join(_full(path), name+'.p')
            saveobj(fullfile, sys.modules[name])


    def reload(self, path):
        """ Load saved objects from disk
        """
        for name in self.names:
            fullfile = join(_full(path), name+'.p')
            sys.modules[name] = loadobj(fullfile)

        self.check()


    def check(self):
        for name in self.names:
            sys.modules[name].check()


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
            raise TypeError("Once defined, parameters cannot be deleted.")
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
        mydict = loadpy(abspath('parameters.py'))
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

    def __setattr__(self, key, val):
        super(SeisflowsPaths, self).__setattr__(*_sub(key, val))

    def load(self):
        mydict = loadpy(abspath('paths.py'))
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
        elif args[1] not in args[0]:
            msg = '%s is not defined.' % args[1]
            super(ParameterError, self).__init__(msg)
        elif key in obj:
            msg = '%s has bad value: ' % args[0], args[1].__getattr__(args[0])
            super(ParameterError, self).__init__(msg)


def custom_import(*names):
    """ Imports SeisFlows module and extracts class of same name. For example,

            custom_import('workflow', 'inversion') 

        imports 'seisflows.workflow.inversion' and, from this module, extracts
        class 'inversion'.
    """
    # parse input arguments
    if len(names) == 0:
        raise Exception(ImportError1)
    if names[0] not in SeisflowsObjects.names:
        raise Exception(ImportError2)
    if len(names) == 1:
        names += (_val(names[0]),)
    if not names[1]:
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
        full_dotted_name = package+'.'+names[0]+'.'+names[1]
        if find_loader(full_dotted_name):
            _exists = True
            break
    if not _exists:
        raise Exception(ImportError3 % 
            (names[0], names[1], names[0].upper()))

    # import module
    module = import_module(full_dotted_name)

    # extract class
    if hasattr(module, names[1]):
        return getattr(module, names[1])
    else:
        raise Exception(ImportError4 % 
            (names[0], names[1], names[1]))


# utility functions

def _path():
    try:
        return SeisflowsPaths()['OUTPUT']
    except:
        cwd = abspath('.')
        return join(cwd, 'output')


def _full(path):
    try:
        return join(abspath(path), 'SeisflowsObjects')
    except:
        raise IOError


def _val(key):
    if key.upper() in SeisflowsParameters():
        return SeisflowsParameters()[key.upper()]
    else:
        return ''


def _sub(key, val):
    if val is None:
        val = ''
    if type(val) not in [str, unicode]:
        raise ValueError('Must be string or unicode: PATH.%s' % key)
    try:
        val = re.sub('~', os.getenv('HOME'), val)
    except:
        raise Exception('Tilde "~" expansion failed.')
    return key, val



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

