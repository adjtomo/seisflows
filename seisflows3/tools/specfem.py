"""
Utilities to interact with, manipulate or call on the external solver, 
i.e., SPECFEM2D/3D/3D_GLOBE
"""
import os
import sys
import numpy as np
import subprocess

from collections import defaultdict
from seisflows3.tools import msg
from seisflows3.tools.math import poissons_ratio
from seisflows3.tools.wrappers import iterable


def call_solver(mpiexec, executable, output="solver.log"):
    """
    Calls MPI solver executable to run solver binaries, used by individual 
    processes to run the solver on system.

    A less complicated version, without error catching, would be
    subprocess.call(f"{mpiexec} {executable}", shell=True)

    :type mpiexec: str
    :param mpiexec: call to mpi. If None (e.g., serial run, defaults to ./)
    :type executable: str
    :param executable: executable function to call
    :type output: str
    :param output: where to redirect stdout
    """
    # mpiexec is None when running in serial mode
    if mpiexec is None:
       exc_cmd = f"./{executable}"

    # Otherwise mpiexec is system dependent (e.g., srun, mpirun)
    else:
        exc_cmd = f"{mpiexec} {executable}"

    try:
        f = open(output, 'w')
        subprocess.check_call(exc_cmd, shell=True, stdout=f)
    except (subprocess.CalledProcessError, OSError):
        print(msg.SolverError.format(exc=exc_cmd))
        sys.exit(-1)
    finally:
        f.close()


class Minmax(defaultdict):
    """
    Keeps track of min, max values of model or kernel
    """
    def __init__(self):
        super(Minmax, self).__init__(lambda: [+np.inf, -np.inf])

    def update(self, keys, vals):
        for key, val in _zip(keys, vals):
            if min(val) < self.dict[key][0]:
                self.dict[key][0] = min(val)
            if max(val) > self.dict[key][1]:
                self.dict[key][1] = max(val)

    def __call__(self, key):
        return self.dict[key]


class Container(defaultdict):
    """
    Dictionary-like object for holding models or kernels
    """
    def __init__(self):
        super(Container, self).__init__(lambda: [])
        self.minmax = Minmax()


def getpar(key, file='DATA/Par_file', sep='=', cast=str):
    """
    Reads parameter from Specfem3D parameter file
    """
    val = None
    with open(file, 'r') as f:
        # Read line by line
        for line in f:
            if line.find(key) == 0:
                # Read key
                key, val = _split(line, sep)
                if not key:
                    continue
                # Read val
                val, _ = _split(val, '#')
                val.strip()
                break

    if val:
        if cast == float:
            val = val.replace('d', 'e')
        return cast(val)
    else:
        print(f"Not found in parameter file: {key}\n")
        raise KeyError


def setpar(key, val, filename="DATA/Par_file", path=".", sep="="):
    """
    Overwrites parameter value to text file. 
    Used to change values in SPECFEM parameter file.

    .. note::
        To ensure we only get the parameter were after, we search for the exact 
        parameter name + first trailing space. Avoids collecting parameters that
        contain other parameter names, e.g., searching for 'SAVE' may return 
        'SAVE_FORWARD', but searching 'SAVE ' should only return the parameter 
        we're after

    .. note::
        Assumes the SPECFEM par file is written in the form: key = value
    """
    val = str(val)

    with open(os.path.join(path, filename), "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines[:]):
            # Search for key + first trailing space
            if line[:len(key) + 1] == f"{key} ":
                _, old_val, *_ = line.strip().split("=")  
                lines[i] = line.replace(old_val, val)
                break

    # Write amended file back to same file name
    with open(os.path.join(path, filename), "w") as f:
        f.writelines(lines)


def check_poissons_ratio(vp, vs, min_val=-1., max_val=0.5):
    """
    Check Poisson's ratio based on Vp and Vs model vectors. Exit SeisFlows3 if
    Poisson's ratio is outside `min_val` or `max_val` which by default are
    set internally by SPECFEM. Otherwise return the value

    :type vp: np.array
    :param vp: P-wave velocity model vector
    :type vs: np.array
    :param vp: S-wave velocity model vector
    :type min_val: float
    :param min_val: minimum model-wide acceptable value for poissons ratio
    :type max_val: float
    :param max_val: maximum model-wide acceptable value for poissons ratio
    :return:
    """
    poissons = poissons_ratio(vp=vp, vs=vs)
    if (poissons.min() < min_val) or (poissons.max() > max_val):
        print(msg.PoissonsRatioError.format(min_val=min_val, max_val=max_val,
                                            pmin=poissons.min(),
                                            pmax=poissons.max())
              )
        sys.exit(-1)
    return poissons


def _split(string, sep):
    """
    Utility function to split a string by a given separation character or str

    :type string: str
    :param string: string to split
    :type sep: str
    :param sep: substring to split by
    """
    n = string.find(sep)
    if n >= 0:
        return string[:n], string[n + len(sep):]
    else:
        return string, ''


def _merge(*parts):
    """
    Utility function to merge various strings together with no breaks
    """
    return ' '.join(parts)


def _zip(keys, vals):
    """
    Zip together keys and vals

    :type keys: dict_keys
    :param keys: keys
    :type vals: dict_values
    :param vals: values
    """
    return zip(iterable(keys), iterable(vals))
