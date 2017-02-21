
from os.path import abspath, join
from shutil import copyfile

import numpy as np


def mread(path, parameters, iproc, prefix='', suffix=''):
    """ Multiparameter read, callable by a single mpi process
    """
    keys = []
    vals = []
    for key in sorted(parameters):
        val = _read1(path, iproc, prefix+key+suffix)
        keys += [key]
        vals += [val]
    return keys, vals


def read(path, parameter, iproc):
    """ Reads from ADIOS container
    """
    raise NotImplementedError


def write(v, path, parameter, iproc):
    """ Writes to ADIOS container
    """
    raise NotImplementedError

