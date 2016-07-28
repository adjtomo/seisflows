
import json
import os
import pickle
import re
import subprocess
import sys
import traceback

from imp import load_source
from importlib import import_module
from os.path import basename

import numpy as np


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
    from seisflows.tools import msg

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


def call_solver_nompi(executable, output='/dev/null'):
    """ Calls non-MPI solver executable

      A less complicated version, without error catching, would be
      subprocess.call(executable, shell=True)
    """
    from seisflows.tools import msg

    try:
        f = open(output,'w')
        subprocess.check_call(
            executable,
            shell=True,
            stdout=f)
    except subprocess.CalledProcessError, err:
        print msg.SolverError % executable
        sys.exit(-1)
    except OSError:
        print msg.SolverError % executable
        sys.exit(-1)
    finally:
        f.close()


def cast(var):
    if isinstance(var, list):
        return var
    elif isinstance(var, float):
        return [var]
    elif isinstance(var, int):
        return [var]
    else:
        raise TypeError


def divides(i, j):
    """Returns true if j divides i"""
    if j is 0:
        return False
    elif i % j:
        return False
    else:
        return True


def exists(name):
    """Wrapper for os.path.exists"""
    if not isinstance(name, basestring):
        return False
    else:
        return os.path.exists(name)


def findpath(name):
    """Resolves absolute path of module"""
    path = import_module(name).__file__

    # adjust file extension
    path = re.sub('.pyc$', '.py', path)

    # strip trailing "__init__.py"
    path = re.sub('__init__.py$', '', path)

    return path


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


def loadpy(abspath):
    # load module
    name = re.sub('.py$', '', basename(abspath))
    module = load_source(name, abspath)

    # strip private attributes
    output = Struct()
    for key, val in vars(module).items():
        if key[0] != '_':
            output[key] = val
    return output


def unique(mylist):
    """Finds unique elements of list"""
    return list(set(mylist))


def loadtxt(filename):
    """Load scalar from text file"""
    return float(np.loadtxt(filename))


def savetxt(filename, v):
    """Save scalar to text file"""
    np.savetxt(filename, [v], '%11.6e')


