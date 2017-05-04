
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
from os.path import basename

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


def call_solver(mpiexec, executable, output='/dev/null'):
    """ Calls MPI solver executable

      A less complicated version, without error catching, would be
      subprocess.call(mpiexec +' '+ executable, shell=True)
    """
    try:
        f = open(output,'w')
        subprocess.check_call(
            mpiexec +' '+ executable,
            shell=True,
            stdout=f)
    except subprocess.CalledProcessError, err:
        print msg.SolverError % (mpiexec +' '+ executable)
        sys.exit(-1)
    except OSError:
        print msg.SolverError % (mpiexec +' '+ executable)
        sys.exit(-1)
    finally:
        f.close()


def divides(i, j):
    """Returns true if j divides i"""
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


