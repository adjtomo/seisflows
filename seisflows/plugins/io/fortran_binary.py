
from os.path import abspath, getsize, join
from shutil import copyfile
from seisflows.tools.tools import iterable

import numpy as np


def read_slice(path, parameters, iproc):
    """ Reads SPECFEM model slice(s)
    """
    vals = []
    for key in iterable(parameters):
        filename = '%s/proc%06d_%s.bin' % (path, iproc, key)
        vals += [_read(filename)]
    return vals


def write_slice(data, path, parameters, iproc):
    """ Writes SPECFEM model slice
    """
    for key in iterable(parameters):
        filename = '%s/proc%06d_%s.bin' % (path, iproc, key)
        _write(data, filename)


def copy_slice(src, dst, iproc, parameter):
    """ Copies SPECFEM model slice
    """
    filename = 'proc%06d_%s.bin' % (iproc, parameter)
    copyfile(join(src, filename), join(dst, filename))


def _read(filename):
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


def _write(v, filename):
    """ Writes Fortran style binary files--data are written as single precision
        floating point numbers
    """
    n = np.array([4*len(v)], dtype='int32')
    v = np.array(v, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        v.tofile(file)
        n.tofile(file)

