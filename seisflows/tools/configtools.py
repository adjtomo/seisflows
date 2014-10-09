
import os as _os
import imp as _imp
import inspect as _inspect
import pickle as _pickle
import sys as _sys
import types as _types

from seisflows.tools import unix
from seisflows.tools.codetools import Struct, abspath, join, loadobj, saveobj


class ParameterObject(object):
    """ Dictionary like object for holding parameters

     Maps the contents of a globally accessible module to a dictionary like
     object that provides some useful utilities as well as optional safegaurds
     against modifying parameters once they are read in.

     Also addresses the following problem: on parallel systems, the environment
     in which jobs are submitted (the head node) may differ from the one in
     which jobs are executed (the compute nodes). To deal with such cases,
     parameters are written to disk by the head node after parameter checking
     so that if necessary they can be read back by the compute nodes.
    """

    def __init__(self,name,path=None):
        self.__dict__['path'] = path
        self.__dict__['locked'] = False
        self.__dict__['module'] = self.load(name)

    def __iter__(self):
        return iter(self.vars.keys())

    def __getattr__(self,key):
        return getattr(self.module,key)

    def __setattr__(self,key,val):
        if  self.__dict__['locked']:
            print 'warning'
            return
        setattr(self.__dict__['module'],key,val)

    def lock(self):
        self.__dict__['locked'] = True

    def unlock(self):
        self.__dict__['locked'] = False

    @property
    def vars(self):
        mystruct = Struct()
        for key,val in vars(self.module).items():
            if key[0] != '_':
                mystruct[key] = val
        return mystruct

    def load(self,name):
        if name in _sys.modules:
            return _sys.modules[name]

        elif self.path:
            # load parameters from disk
            mydict = loadobj(join(self.path,name+'.p'))

            # register module
            module = _types.ModuleType(name)
            _sys.modules[name] = module

            # populate module
            for key,val in mydict.items():
                if key[0] != '_':
                    setattr(module,key, val)
            return module

        else:
            raise Exception


def getclass(*args):
    """ Given name of module realtive to package directory, returns
     corresponding class
    """
    if not args[-1]:
        return object # return dummy class

    # first, try importing relative to main package directory
    list = _parse(args,package='seisflows')
    string = '.'.join(list)
    if _exists(list):
        module = getattr(_import(string),list[-1])
        return module

    # next, try importing relative to extensions directory
    list = _parse(args,package='seisflows.extensions')
    string = '.'.join(list)
    if _exists(list):
        module = getattr(_import(string),list[-1])
        return module

    # last, try importing relative to working directory
    list = _parse(args)
    string = '.'.join(list)
    if _exists(list):
        module = getattr(_import(string),list[-1])
        return module

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
        try:
            string = '.'.join(_parse([obj],package='seisflows'))
            module = _import(string)
            path = module.__file__
        except:
            raise Exception
    return _os.path.dirname(path)



### utility functions

def _import(string):
    "Imports module from string"
    return __import__(string,fromlist='dummy')


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
