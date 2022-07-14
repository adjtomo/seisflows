"""
General utility functions that are mostly concerend with file manipulation,
but also math and calling functions as well.
"""
import os
import re
import yaml
import numpy as np
from seisflows.core import Dict
from seisflows import logger


def log_status(func):
    """
    Decorator function that logs the completion status of a function to a
    state file. This is used for checkpointing a workflow and resuming
    failed workflows without repeating computational intense tasks
    """
    STATE_FILE = os.path.join(os.getcwd(), "sfstatefile")

    def logged_func():
        """Log the completion status of the function"""
        try:
            func()
            with open(STATE_FILE, "a") as f:
                f.write(f"{func.__name__}\tCOMPLETED")
        except Exception as e:
            f.write(f"{func.__name__}\tFAILED")
            logger.error(e)
            raise

    lines = open(STATE_FILE, "r").readlines()
    for line in lines:
        function, status = line.split(" ")
        if func.__name__ == function:
            if status == "COMPLETE":
                return
            elif status == "FAILED":
                return logged_func()
        else:
            return logged_func()


def set_task_id(task_id):
    """
    Set the SEISFLOWS_TASKID in os environs

    .. note::
        Mostly used for debugging/testing purposes as a way of mimicing
        system.run() assigning task ids to child processes

    :type task_id: int
    :param task_id: integer task id to assign to the current working environment
    """
    os.environ["SEISFLOWS_TASKID"] = str(task_id)


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

