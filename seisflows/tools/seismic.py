"""
Seismic signal processing related tools
"""
import os
import sys
import numpy as np
import subprocess

from collections import defaultdict
from seisflows.tools import msg, unix
from seisflows.tools.tools import iterable


def call_solver(mpiexec, executable, output='solver.log'):
    """
    Calls MPI solver executable

    A less complicated version, without error catching, would be
    subprocess.call(f"{mpiexec} {executable}", shell=True)

    :type mpiexec: str
    :param mpiexec: call to mpi
    :type executable: str
    :param executable: executable function to call
    :type output: str
    :param output: where to redirect stdout
    """
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


class Writer(object):
    """
    Utility for appending values to text files.
    Used for writing statistical outputs
    """
    def __init__(self, path="./output.stat"):
        self.path = os.path.abspath(path)
        try:
            os.mkdir(path)
        except FileExistsError:
            print(f"Warning, {os.path.basename(path)} exists\n"
                  "Appending to this files without deleting them may lead to "
                  "unintended consequences")

        self.__call__('step_count', 0)

    def __call__(self, filename, val):
        """
        Define function for calling Writer

        :type filename: str
        :param filename: file to write to
        :type val: float
        :param val: value to write to file
        """
        fullfile = os.path.join(self.path, filename)

        with open(fullfile, 'a') as f:
            f.write(f"{val:e}\n")


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


def setpar(key, val, filename='DATA/Par_file', path='.', sep='='):
    """
    Writes parameter to text file. Used to change Specfem3D parameter file
    """
    val = str(val)

    # Read line by line
    with open(os.path.join(path, filename), "r") as file:
        lines = []
        for line in file:
            if line.find(key) == 0:
                # Read key
                key, _ = _split(line, sep)
                # Read comment
                _, comment = _split(line, '#')
                n = len(line) - len(key) - len(val) - len(comment) - 2
                # Replace line with new value
                if comment:
                    line = _merge(key, sep, val, ' '*n, '#', comment)
                else:
                    line = _merge(key, sep, str(val), '\n')
            lines.append(line)

    # Write file
    with open(os.path.join(path, filename), 'w') as file:
        file.writelines(lines)


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
    return ''.join(parts)


def _zip(keys, vals):
    """
    Zip together keys and vals

    :type keys: dict_keys
    :param keys: keys
    :type vals: dict_values
    :param vals: values
    """
    return zip(iterable(keys), iterable(vals))
