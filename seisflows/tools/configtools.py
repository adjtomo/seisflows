
import os as _os
import imp as _imp
import inspect as _inspect
import pickle as _pickle
import sys as _sys
import types as _types

from seisflows.tools import unix
from seisflows.tools.codetools import abspath, exists ,join, loadjson, loadobj

class Config:
    pass


class GlobalStruct(object):
    """ Globally accessible object for holding parameters
    """

    def __new__(self,name,path='.'):
        if name in _sys.modules:
            return _sys.modules[name]
        else:
            return object.__new__(self)

    def __init__(self,name,path='.'):
        if name not in _sys.modules:
            self.fromfile(name,path)
            _sys.modules[name] = self

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

    def fromfile(self,name,path):
        try:
            path = abspath(path)
        except:
            raise IOError, path

        fullname = join(path,name)
        if exists(fullname+'.p'):
            vars = loadobj(fullname+'.p')
        elif exists(fullname+'.json'):
            vars = loadjson(fullname+'.json')
        elif exists(fullname+'.py'):
            vars = _vars(_import(name,path))
        else:
            raise IOError, fullname

        super(GlobalStruct,self).__setattr__('__dict__',vars)




def getclass(*args):
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

    # last, try importing relative to working directory
    list = _parse(args)
    string = '.'.join(list)
    if _exists(list):
        moduleobj = getattr(_import(string),list[-1])
        return moduleobj

    raise Exception


def getmodule(obj):
    """Tries to determine module in which class, function, or method is defined
   """
    if isinstance(obj,_types.FunctionType):
        name = obj.__module__

    elif isinstance(obj,_types.MethodType):
        name = obj.im_class.__module__

    elif isinstance(obj,_types.InstanceType):
        name = obj.__class__.__module__

    else:
        name = _inspect.getmodule(obj).__name__

    if name == '__main__':
        name,_ = _os.path.splitext(
                   _os.path.basename(
                     _inspect.getfile(obj)))
    return name


def getpath(obj):
    """ Determines absolute path of module from name, relative path, or
      ModuleObject
    """
    if isinstance(obj,_types.ModuleType):
        path = obj.__file__
    elif _os.path.isfile(obj):
        path = _os.path.abspath(obj)
    else:
        string = '.'.join(_parse([obj],package='seisflows'))
        moduleobj = _import(string)
        path = moduleobj.__file__
    return _os.path.dirname(path)



### utility functions

def _import(string,path=None):
    "Imports module from string"
    if path:
        # temporarily adjust python path
        _sys.path.append(path)

    moduleobj = __import__(string,fromlist='dummy')

    if path:
        _sys.path.pop()

    return moduleobj


def _vars(obj):
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
    "Check if a module exists without importing it"
    try:
        path = None
        for part in parts[:-1]:
            args = _imp.find_module(part,path)
            obj = _imp.load_module(part,*args)
            path = obj.__path__
        args = _imp.find_module(parts[-1],path)
        return True
    except ImportError:
        return False

