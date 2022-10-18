#!/usr/bin/env python3
"""
Only required when system==cluster (or any subclass of cluster)

This script is a wrapper for running tasks on systems during an active workflow.
Loads in functions/methods and their keyword arguments from pickle files which
have been written by system.run(). Runs functions in order of list.

.. note::
    Not to be called by the user, this script is to be called by system.run()

.. rubric::
    >> python run --funcs function_list.p --kwargs kwarg_dict.p
"""
import argparse
import dill
import os

from seisflows import logger
from seisflows.tools.config import config_logger


def parse_args():
    """
    Get command line arguments
    """
    parser = argparse.ArgumentParser("Run arguments for system submitted tasks")

    parser.add_argument("-f", "--funcs", type=str, nargs="?", required=True,
                        help="path to pickle file containing a list of "
                             "functions/methods that should be run by the "
                             "submitted process"
                        )
    parser.add_argument("-k", "--kwargs", type=str, nargs="?", required=False,
                        default=None,
                        help="path to pickle file containing a dictionary of "
                             "keyword argumnets that should be passed to the "
                             "functions")
    parser.add_argument("-e", "--environment", type=str, nargs="?",
                        required=False,
                        help="Optional comma-separated environment variables, "
                             "which should be given as "
                             "VARNAME1=value1,VARNAME2=value2 and so on. These "
                             "will be separated and instantiated into Python's "
                             "os.environ")

    return parser.parse_args()


def export(myenv):
    """
    Exports comma delimited list of environment variables also allows deleting 
    environment variables by providing VARNAME with no corresponding value 

    e.g. VARNAME1=value1,VARNAME2=value2,VARNAME3
    will add VARNAME1 and VARNAME2 to the environment with corresponding values, 
    and remove VARNAME3 from the environment

    .. note::
        The ability to delete environment variables came from the Maui upgrade
        to Slurm 21.08, which enforced mutually exclusivity of --mem-per-cpu 
        and --mem-per-node, which are both defined on cross-cluster submissions.
        We needed a mechanism to remove one of these

    :type myenv: str
    :param myenv: the system environment to take variables from
    """
    for item in myenv.split(","):
        if item:
            try:
                key, val = item.split("=")
                os.environ[key] = val
            # Variables to be deleted will not split on '=', throwing ValueError
            except ValueError:
                try:
                    del os.environ[item]
                # If a NoneType sneaks through, it will throw TypeEror on 'del'
                except TypeError:
                    continue


if __name__ == '__main__':
    """Runs task within a currently executing workflow """
    args = parse_args()

    if args.environment:
        export(args.environment)

    # Load the functions
    with open(args.funcs, "rb") as f:
        funcs = dill.load(f)

    # Load the kwargs. Optional, if No kwargs then funcs will be run bare
    if args.kwargs is not None:
        with open(args.kwargs, "rb") as f:
            kwargs = dill.load(f)
    else:
        kwargs = {}

    # Hold onto the main handler for after the task is finished
    handlers = logger.handlers

    # Redirect the logger to point at stdout, which the Cluster system (or sub
    # class of it) will redirect to a task specific log file. Currently no
    # way of accessing user defiend 'verbose' or 'log_level' parameters
    config_logger(filename=None, stream_handler=True)

    # Evaluate the function with given keyword arguments
    for func in funcs:
        func(**kwargs)

    # Replace original log handler now that the tasks have finished
    logger.handlers = handlers


