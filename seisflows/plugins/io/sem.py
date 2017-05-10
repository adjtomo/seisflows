
from os.path import abspath, getsize, join
from shutil import copyfile
from seisflows.tools.tools import iterable

import numpy as np


def read(path, parameters, iproc):
    """ Reads SPECFEM database file(s)
    """
    vals = []
    for key in iterable(parameters):
        filename = '%s/proc%06d_%s.bin' % (path, iproc, key)
        vals += [_read_bin(filename)]
    return vals


def write(data, path, parameter, iproc):
    """ Writes a single SPECFEM database file
    """
    filename = 'proc%06d_%s.bin' % (iproc, parameter)
    _write_bin(data, join(path, filename))


def copy(src, dst, iproc, parameter):
    """ Copies SPECFEM database file
    """
    filename = 'proc%06d_%s.bin' % (iproc, parameter)
    copyfile(join(src, filename), join(dst, filename))


def _read_bin(filename):
    """ Reads Fortran style binary data into numpy array
    """
    nbytes = getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n==nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype='float32')
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype='float32')
            return data


def _write_bin(v, filename):
    """ Writes Fortran style binary files--data are written as single precision
        floating point numbers
    """
    n = np.array([4*len(v)], dtype='int32')
    v = np.array(v, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        v.tofile(file)
        n.tofile(file)

