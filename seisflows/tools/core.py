"""
Core utility functions and classes which help define the working structure of
SeisFlows pacakge.
"""
import os
import re
import yaml
import numpy as np
from seisflows import logger

# Acceptable environment variables assigned to individually running tasks when
# running SeisFlows on a system
ENV_VARIABLES = ["SEISFLOWS_TASKID", "SLURM_ARRAY_TASK_ID"]


class Dict(dict):
    """
    A dictionary replacement which allows for easier parameter access through
    getting and setting attributes. Also has some functionality to make string
    printing prettier
    """
    def __str__(self):
        """Pretty print dictionaries and first level nested dictionaries"""
        str_ = ""
        try:
            longest_key = max([len(_) for _ in self.keys()])
            for key, val in self.items():
                str_ += f"{key:<{longest_key}}: {val}\n"
        except ValueError:
            pass
        return str_

    def __repr__(self):
        """Pretty print when calling an instance of this object"""
        return self.__str__()

    def __getattr__(self, key):
        """Attribute-like access of the internal dictionary attributes"""
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"{key} not found in Dict")

    def __setattr__(self, key, val):
        """Setting attributes can only be performed one time"""
        self.__dict__[key] = val


class Null:
    """
    A null object that always and reliably does nothing
    """
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        return self

    def __bool__(self):
        return False

    def __nonzero__(self):
        return False

    def __getattr__(self, key):
        return self

    def __setattr__(self, key, val):
        return self

    def __delattr__(self, key):
        return self


def get_task_id():
    """
    Task IDs are assigned to each child process spawned by the system module
    during a SeisFlows workflow. SeisFlows modules use this Task ID to keep
    track of embarassingly parallel process, e.g., solver uses the Task ID to
    determine which source is being considered.

    :rtype: int
    :return: task id for given solver
    """
    for env_var in ENV_VARIABLES:
        _taskid = os.getenv(env_var)
        if _taskid is not None:
            return int(_taskid)
    else:
        logger.warning("Environment Task ID variable not found. Assigning 0")
        return 0


def set_task_id(task_id):
    """
    Set the SEISFLOWS_TASKID in os environs for local workflows. If running
    on HPC systems, running array jobs will assign the Task ID

    .. note::
        Mostly used for debugging/testing purposes as a way of mimicing
        system.run() assigning task ids to child processes

    :type task_id: int
    :param task_id: integer task id to assign to the current working environment
    """
    os.environ["SEISFLOWS_TASKID"] = str(task_id)


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
        Replace statement is catch-all, so we assume that there is only one
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

