"""
Functions to read and write FORTRAN binary files that are outputted by Specfem
"""
import os
import numpy as np
from shutil import copyfile
from seisflows.tools.utils import iterable


def read_slice(path, parameters, iproc):
    """ 
    Reads SPECFEM model slice(s)
    
    :type path: str
    :param path: path to the database files
    :type parameters: str
    :param parameters: parameters to read, e.g. 'vs', 'vp'
    :type iproc: int
    :param iproc: processor/slice number to read
    """
    vals = []
    for key in iterable(parameters):
        filename = os.path.join(path, f"proc{int(iproc):06d}_{key}.bin")
        vals += [_read(filename)]
    return vals


def write_slice(data, path, parameters, iproc):
    """ 
    Writes SPECFEM model slice

    :type data: seisflows.Container
    :param data: data to be written to a slice
    :type path: str
    :param path: path to the database files
    :type parameters: str
    :param parameters: parameters to write, e.g. 'vs', 'vp'
    :type iproc: int
    :param iproc: processor/slice number to write
    """
    for key in iterable(parameters):
        filename = os.path.join(path, f"proc{int(iproc):06d}_{key}.bin")
        _write(data, filename)


def copy_slice(src, dst, iproc, parameter):
    """ 
    Copies SPECFEM model slice

    :type src: str
    :param src: source location to copy slice from
    :type dst: str
    :param dst: destination location to copy slice to
    :type parameter: str
    :param parameter: parameters to copy, e.g. 'vs', 'vp'
    :type iproc: int
    :param iproc: processor/slice number to copy
    """
    filename = f"proc{int(iproc):06d}_{parameter}.bin"
    copyfile(os.path.join(src, filename), 
             os.path.join(dst, filename))


def _read(filename):
    """ 
    Reads Fortran style binary data into numpy array
    """
    nbytes = os.path.getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n == nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype='float32')
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype='float32')
            return data


def _write(v, filename):
    """ 
    Writes Fortran style binary files
    Data are written as single precision floating point numbers

    .. note::
        FORTRAN unformatted binaries are bounded by an INT*4 byte count. This
        function mimics that behavior by tacking on the boundary data.
        https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc4/index.html
    """
    n = np.array([4 * len(v)], dtype='int32')
    v = np.array(v, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        v.tofile(file)
        n.tofile(file)

