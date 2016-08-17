
from os.path import abspath, getsize, join
from shutil import copyfile

import numpy as np


def mread(path, parameters, iproc, prefix='', suffix=''):
    """ Multiparameter read, callable by a single mpi process
    """
    keys = []
    vals = []
    for key in sorted(parameters):
        val = read(path, iproc, prefix+key+suffix)
        keys += [key]
        vals += [val]
    return keys, vals


def read(path, parameter, iproc):
    """ Reads a single SPECFEM database file
    """
    filename = 'proc%06d_%s.bin' % (iproc, parameter)
    return _read_bin(join(path, filename))


def write(v, path, parameter, iproc):
    """ Writes a single SPECFEM database file
    """
    filename = 'proc%06d_%s.bin' % (iproc, parameter)
    _write_bin(v, join(path, filename))


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
            v = np.fromfile(file, dtype='float32')
            return v[:-1]
        else:
            file.seek(0)
            v = np.fromfile(file, dtype='float32')
            return v


def _write_bin(v, filename):
    """ Writes Fortran style binary data--data are written as single precision
        floating point numbers
    """
    n = np.array([4*len(v)], dtype='int32')
    v = np.array(v, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        v.tofile(file)
        n.tofile(file)

