"""
Functions to read and write ASCII model (.dat) files used by SPECFEM2D
"""
import os
import numpy as np
from glob import glob
from shutil import copyfile


def read_slice(path, parameters, iproc):
    """ 
    Reads SPECFEM model slice(s) based on .dat ASCII files
    
    :type path: str
    :param path: path to the database files
    :type parameters: str
    :param parameters: parameters to read, e.g. 'vs', 'vp'
    :type iproc: int
    :param iproc: processor/slice number to read
    :rtype: list of np.array
    :return: list of arrays corresponding to model parameters in given order
    """
    filename = _get_filename(path, iproc)
    available_parameters = _get_available_parameters(filename)
    model = np.loadtxt(filename).T
    
    vals = []
    for key in parameters:
        vals += [model[available_parameters.index(key)]]
    
    return vals


def write_slice(data, path, parameters, iproc):
    """ 
    Writes SPECFEM model slice

    !!! This won't work because we need access to the spatial components that
    !!! are only the model

    :type data: seisflows.Container
    :param data: data to be written to a slice
    :type path: str
    :param path: path to the database files
    :type parameters: str
    :param parameters: parameters to write, e.g. 'vs', 'vp'
    :type iproc: int
    :param iproc: processor/slice number to write
    """
    for key in parameters:
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
    filename = os.path.basename(_get_filename(src, iproc))
    copyfile(os.path.join(src, filename), 
             os.path.join(dst, filename))


def _get_filename(path, iproc):
    """
    ASCII .dat files list the available parameters in the fileid, meaning
    there is no standard format for retrieving files. Use glob to search for
    the file based on file extension.

    :type path: str
    :param path: path to the database files
    :type iproc: int
    :param iproc: processor/slice number to read
    :rtype: str
    :return: filename of the model
    """
    filename_glob = os.path.join(path, f"proc{int(iproc):06d}_*.dat")
    filename = glob(filename_glob)
    assert(len(filename) == 1), \
        f"Expected only one .dat file, found {len(filename)}"

    return filename[0]


def _get_available_parameters(filename):
    """
    The available parameters are listed in the file name. Split off the 
    uncessary text and return the listend parameters.

    :type filename: str
    :param filename: filename to check parameters from
    :rtype: list
    :return: list of parameters from the file id
    """
    fid = os.path.splitext(os.path.basename(filename))[0]
    _, *available_parameters = fid.split("_")
    
    return available_parameters


def _write(v, filename):
    """ 
    Writes Fortran style binary files
    Data are written as single precision floating point numbers
    """
    n = np.array([4 * len(v)], dtype='int32')
    v = np.array(v, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        v.tofile(file)
        n.tofile(file)

