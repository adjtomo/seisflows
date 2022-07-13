"""
General utility functions that are mostly concerend with file manipulation,
but also math and calling functions as well.
"""
import os
import re
import time
import yaml
import subprocess
import numpy as np
from importlib import import_module
from pkgutil import find_loader
from seisflows.core import Dict
from seisflows import logger


def get_task_id():
    """
    Task IDs are assigned to each child process spawned by the system module
    during a SeisFlows workflow. SeisFlows modules use this Task ID to keep
    track of embarassingly parallel process, e.g., solver uses the Task ID to
    determine which source is being considered.

    :rtype: int
    :return: task id for given solver
    """
    _taskid = os.getenv("SEISFLOWS_TASKID")
    if _taskid is None:
        _taskid = 0
        logger.warning("Environment variable 'SEISFLOWS_TASKID' not found. "
                       "Assigning Task ID == 0")
    return int(_taskid)


def load_yaml(filename):
    """
    Define how the PyYaml yaml loading function behaves.
    Replaces None and inf strings with NoneType and numpy.inf respectively

    :type filename: str
    :param filename: .yaml file to load in
    :rtype: Dict
    :return: Dictionary containing all parameters in a YAML file
    """
    # work around PyYAML bugs
    yaml.SafeLoader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.'))

    with open(filename, 'r') as f:
        mydict = Dict(yaml.safe_load(f))

    if mydict is None:
        mydict = Dict()

    # Replace 'None' and 'inf' values to match expectations
    for key, val in mydict.items():
        if val == "None":
            mydict[key] = None
        if val == "inf":
            mydict[key] = np.inf

    return mydict


def iterable(arg):
    """
    Make an argument iterable

    :param arg: an argument to make iterable
    :type: list
    :return: iterable argument
    """
    if not isinstance(arg, (list, tuple)):
        return [arg]
    else:
        return arg


def number_fid(fid, i=0):
    """
    Number a filename. Used to store old log files without overwriting them.
    Premise is, if you have a file e.g., called: output.txt
    This function would return incrementing filenames:
    output_000.txt, output_001.txt, output_002.txt, ouput_003.txt ...

    .. note::
        Replace statement is catch all so we assume that there is only one \
        instance of the file extension in the entire path.

    :type fid: str
    :param fid: path to the file that you want to increment
    :type i: int
    :param i: number to append to file id
    :rtype: str
    :return: filename with appended number. filename ONLY, will strip away
        the original path location
    """
    fid_only = os.path.basename(fid)
    ext = os.path.splitext(fid_only)[-1]  # e.g., .txt
    new_ext = f"_{i:0>3}{ext}"   # e.g., _000.txt
    new_fid = fid_only.replace(ext, new_ext)
    return new_fid


def nproc():
    """
    Get the number of processors available

    :rtype: int
    :return: number of processors
    """
    try:
        return _nproc_method1()
    except EnvironmentError:
        return _nproc_method2()


def _nproc_method1():
    """
    Used subprocess to determine the number of processeors available

    :rtype: int
    :return: number of processors
    """
    # Check if the command `nproc` works
    if not subprocess.getstatusoutput('nproc')[0] == 0:
        raise EnvironmentError

    num_proc = int(subprocess.getstatusoutput('nproc')[1])

    return num_proc


def _nproc_method2():
    """
    Get number of processors using /proc/cpuinfo

    Bryant: This doesnt work?

    :rtype: int
    :return: number of processors
    """
    if not os.path.exists('/proc/cpuinfo'):
        raise EnvironmentError

    stdout = subprocess.check_output(
        "cat /proc/cpuinfo | awk '/^processor/{print $3}'", shell=True)
    num_proc = len(stdout.split('\n'))

    return num_proc

