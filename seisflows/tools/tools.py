
import json
import os
import pickle
import re
import subprocess
import sys
import time
import traceback

from imp import load_source
from importlib import import_module
from pkgutil import find_loader
from os.path import basename, exists
from subprocess import check_output

import numpy as np

from seisflows.tools import msg


class Struct(dict):
    def __init__(self, *args, **kwargs):
        super(Struct, self).__init__(*args, **kwargs)
        self.__dict__ = self


def call(*args, **kwargs):
    if 'shell' not in kwargs:
        kwargs['shell'] = True
    subprocess.check_call(*args, **kwargs)


def diff(list1, list2):
    """ Difference between two lists
    """
    c = set(list1).union(set(list2))
    d = set(list1).intersection(set(list2))
    return list(c - d)


def divides(i, j):
    """True if j divides i"""
    if j is 0:
        return False
    elif i % j:
        return False
    else:
        return True


def exists(names):
    """Wrapper for os.path.exists"""
    for name in iterable(names):
        if not name:
            return False
        elif not isinstance(name, basestring):
            raise TypeError
        elif not os.path.exists(name):
            return False
    else:
        return True


def findpath(name):
    """Resolves absolute path of module"""
    path = import_module(name).__file__

    # adjust file extension
    path = re.sub('.pyc$', '.py', path)

    # strip trailing "__init__.py"
    path = re.sub('__init__.py$', '', path)

    return path


def iterable(arg):
    if not isinstance(arg, (list, tuple)):
        return [arg]
    else:
        return arg


def module_exists(name):
    return find_loader(name)


def package_exists(name):
    return find_loader(name)


def pkgpath(name):
    for path in import_module('seisflows').__path__:
        if name+'/seisflows' in path:
            return path


def timestamp():
    return time.strftime('%H:%M:%S')


def loadobj(filename):
    """Load object using pickle"""
    with open(filename, 'rb') as file:
        return pickle.load(file)


def saveobj(filename, obj):
    """Save object using pickle"""
    with open(filename, 'wb') as file:
        pickle.dump(obj, file)


def loadjson(filename):
    """Load object using json"""
    with open(filename, 'rb') as file:
        return json.load(file)


def savejson(filename, obj):
    """Save object using json"""
    with open(filename, 'wb') as file:
        json.dump(obj, file, sort_keys=True, indent=4)


def loadpy(filename):
    if not exists(filename):
        print msg.FileError % filename
        raise IOError

    # load module
    name = re.sub('.py$', '', basename(filename))
    module = load_source(name, filename)

    # strip private attributes
    output = Struct()
    for key, val in vars(module).items():
        if key[0] != '_':
            output[key] = val
    return output


def loadnpy(filename):
    """Loads numpy binary file."""
    return np.load(filename)


def savenpy(filename, v):
    """Saves numpy binary file."""
    np.save(filename, v)
    os.rename(filename + '.npy', filename)



def loadyaml(filename):
    import yaml

    with open(filename, 'rb') as file:
        dict = yaml.load(file)

    # replace None
    if 'None' in dict.values():
        for key,val in dict.items():
            if val=='None': dict[key]=None

    return dict


def getset(arg):
    if not arg:
        return set()
    elif isinstance(arg, basestring):
        return set([arg])
    else:
        return set(arg)


def loadtxt(filename):
    """Load scalar from text file"""
    return float(np.loadtxt(filename))


def savetxt(filename, v):
    """Save scalar to text file"""
    np.savetxt(filename, [v], '%11.6e')


def nproc():
    try:
        return _nproc1()
    except:
        return _nproc2()


def _nproc1():
    # get number of processors using nproc
    if not which('nproc'):
        raise EnvironmentError
    stdout = check_output('nproc --all', shell=True)
    nproc = int(stdout.strip())
    return nproc


def _nproc2():
    # get number of processors using /proc/cpuinfo
    if not exists('/proc/cpuinfo'):
        raise EnvironmentError
    stdout = check_output("cat /proc/cpuinfo | awk '/^processor/{print $3}'", 
                shell=True)
    nproc = len(stdout.split('\n'))
    return nproc

