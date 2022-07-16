"""
General utility functions that are mostly concerend with file manipulation,
but also math and calling functions as well.
"""
import os
import re
import time
import yaml
import numpy as np
from seisflows.core import Dict
from seisflows import logger


class TaskIDError(Exception):
    """
    A specific error that gets called when tasks are not run on system,
    i.e., when we can't find 'SEISFLOWS_TASKID' in the environment variables.
    This means we are attempting to access child process variables inside
    the parent process.
    """
    pass


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


def get_task_id(force=False):
    """
    Task IDs are assigned to each child process spawned by the system module
    during a SeisFlows workflow. SeisFlows modules use this Task ID to keep
    track of embarassingly parallel process, e.g., solver uses the Task ID to
    determine which source is being considered.

    :type force: bool
    :param force: If no task id is found, force set it to 0
    :rtype: int
    :return: task id for given solver
    :raises TaskIDError: if no environment variable is found
    """
    _taskid = os.getenv("SEISFLOWS_TASKID")
    if _taskid is None:
        if force:
            _taskid = 0
            logger.warning("Environment variable 'SEISFLOWS_TASKID' not found. "
                           "Assigning Task ID == 0")
        else:
            raise TaskIDError("Environment variable 'SEISFLOWS_TASKID' not "
                              "found. Please make sure the process asking "
                              "for task id is called by system.")
    return int(_taskid)


def log_status(func):
    """
    Decorator function that logs the completion status of a function to a
    state file. This is used for checkpointing a workflow and resuming
    failed workflows without repeating computational intense tasks

    :type func: function
    """
    raise NotImplementedError("This is not working as expected")

    STATE_FILE = os.path.join(os.getcwd(), "sfstatefile")

    if not os.path.exists(STATE_FILE):
        with open(STATE_FILE, "w") as f:
            f.write(f"# SeisFlows State File\n")
            f.write(f"# {time.asctime()}\n")
            f.write(f"# =========================\n")

    def logged_func():
        """Log the completion status of the function"""
        try:
            output = func()
            with open(STATE_FILE, "a") as f:
                f.write(f"{func.__name__}\tCOMPLETED\n")
        except Exception as e:
            with open(STATE_FILE, "a") as f:
                f.write(f"{func.__name__}\tFAILED\n")
            logger.error(e)
            raise
        return output

    lines = open(STATE_FILE, "r").readlines()
    for line in lines:
        if line.startswith("#"):
            continue
        function, status = line.strip().split("\t")
        if func.__name__ == function:
            if status == "COMPLETED":
                return
            elif status == "FAILED":
                return logged_func()
        else:
            return logged_func()
    else:
        return logged_func()


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

