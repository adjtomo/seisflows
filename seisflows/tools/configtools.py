
import os
import imp
import inspect
import sys
import types

from seisflows.tools import unix
from seisflows.tools.codetools import Struct, abspath, exists, glob, join, loadjson, loadobj, savejson, saveobj


class ConfigObj(object):
    """ Makes objects globally accessible by registering them in sys.modules

        Also provides methods for reading and writing objects to disk
    """

    def __init__(self,name):
        if name not in sys.modules:
             sys.modules[name] = set()

        self.name = name
        self.keys = sys.modules[name]

    def register(self,key,val):
        "Registers an object"
        sys.modules[key] = val
        self.keys.add(key)

    def unregister(self,key):
        "Unregisters an object"
        sys.modules.__delattr__(key)
        self.keys.remove(key)

    def save(self,name,path='.'):
        "Saves current state"
        try:
            fullpath = join(abspath(path),name)
        except:
            raise IOError, path

        unix.mkdir(fullpath)
        for key in self.keys:
            saveobj(fullpath+'/'+key+'.p',sys.modules[key])

    def load(self,name,path='.'):
        try:
            fullpath = join(abspath(path),name)
        except:
            raise IOError, path

        files = glob(join(fullpath,'*.p'))
        for file in files:
            key,_ = os.path.splitext(unix.basename(file))
            self.keys.add(key)
            sys.modules[key] = loadobj(file)



class ParameterObj(object):
    """ Dictionary like object for holding parameters

        Makes parameters globally accessible by registering itself in sys.modules
    """

    def __new__(self,name,path='.'):
        if name in sys.modules:
            return sys.modules[name]
        else:
            return object.__new__(self)

    def __init__(self,name,path='.'):
        if name not in sys.modules:
            sys.modules[name] = self

    def __iter__(self):
        return iter(self.__dict__.keys())

    def __getattr__(self,key):
        return self.__dict__[key]

    def __setattr__(self,key,val):
        if key in self.__dict__:
            raise TypeError
        self.__dict__[key] = val

    def __delattr__(self,key):
        if key in self.__dict__:
            raise TypeError
        raise KeyError

    def update(self,newdict):
        super(ParameterObj,self).__setattr__('__dict__',newdict)

    def save(self,name):
        savejson(name,self.__dict__)


def loadclass(*args):
    """ Given name of module relative to package directory, returns
     corresponding class object
    """
    if not args[-1]:
        return object # return dummy class

    # first, try importing relative to main package directory
    list = _parse(args,package='seisflows')
    string = '.'.join(list)
    if _exists(list):
        moduleobj = getattr(_import(string),list[-1])
        return moduleobj

    # next, try importing relative to extensions directory
    list = _parse(args,package='seisflows.extensions')
    string = '.'.join(list)
    if _exists(list):
        moduleobj = getattr(_import(string),list[-1])
        return moduleobj

    raise ImportError


def loadvars(*args,**kwargs):
    return _vars(_import(*args,**kwargs))



def findpath(obj):
    """ Determines absolute path of module from name, relative path, or
      ModuleObject
    """
    if isinstance(obj,types.ModuleType):
        path = obj.__file__
    elif os.path.isfile(obj):
        path = os.path.abspath(obj)
    else:
        string = '.'.join(_parse([obj],package='seisflows'))
        moduleobj = _import(string)
        path = moduleobj.__file__
    return os.path.dirname(path)



### utility functions

def _import(string,path=None):
    "Imports from string"
    if path:
        # temporarily adjust python path
        sys.path.append(path)

    moduleobj = __import__(string,fromlist='dummy')

    if path:
        sys.path.pop()

    return moduleobj


def _vars(obj):
    "Returns an object's __dict__ with private varibles removed"
    mydict = {}
    for key,val in vars(obj).items():
        if key[0] != '_':
            mydict[key] = val
    return mydict


def _parse(args,package=None):
    arglist = []
    for arg in args:
        arglist.extend(arg.split('.'))
    if package:
        parglist = package.split('.')
        nn = min(len(arglist),len(parglist))
        for ii in range(nn+1):
            if parglist[ii:nn] == arglist[0:nn-ii]:
                break
        arglist = parglist[0:ii] + arglist
    return arglist


def _exists(parts):
    "Checks if a module exists without importing it"
    try:
        path = None
        for part in parts[:-1]:
            args = imp.find_module(part,path)
            obj = imp.load_module(part,*args)
            path = obj.__path__
        args = imp.find_module(parts[-1],path)
        return True
    except ImportError:
        return False

