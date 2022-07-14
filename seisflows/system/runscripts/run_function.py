#!/usr/bin/env python3
"""
Only required when system==cluster (or any subclass of cluster)

This script is a wrapper for running tasks on systems during an active workflow.
Acts as a Python script to submit certain SeisFlows functions or tasks to a
compute system.

.. note::
    Not to be called by the user, this script is to be called by system.run()

.. rubric::
    >> python run --output ./OUTPUT --classname solver \
        --funcname eval_func
    OR
    >> sbatch run --output ./OUTPUT --classname solver --funcname eval_func
"""
import os
import sys
import pickle
import argparse

from seisflows.config import load, config_logger


def parse_args():
    """
    Get command line arguments
    """
    parser = argparse.ArgumentParser("Run arguments for system submitted tasks")
    parser.add_argument("-o", "--output", type=str, nargs="?", required=True,
                        help="the SeisFlows output directory used to load the "
                             "active working state from inside the compute node"
                        )
    parser.add_argument("-c", "--classname", type=str, nargs="?", required=True,
                        help="the SeisFlows class from within which the "
                             "desired function is defined. Available options "
                             "are defined in seisflows.config.NAMES"
                        )
    parser.add_argument("-f", "--funcname", type=str, nargs="?", required=True,
                        help="the function name from the chosen `classname`. "
                             "This function will be executed on the compute "
                             "node")
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
        try:
            key, val = item.split("=") 
            os.environ[key] = val
        # Variables to be deleted will not split on '=', throwing the ValueError
        except ValueError:
            del os.environ[item]


if __name__ == '__main__':
    """ 
    Runs task within a currently executing workflow 
    """
    args = parse_args()

    if args.environment:
        export(args.environment)

    # Load the last checkpointed working state from the 'seisflows_?.p` files
    # Allowing access through sys.modules
    load(args.output)

    # Load keyword arguments required by this function
    # Files will be something like: 'solver_eval_func.p'
    kwargs_fid = f"{args.classname}_{args.funcname}.p"
    kwargs_path = os.path.join(args.output, "kwargs", kwargs_fid)
    with open(kwargs_path, "rb") as f:
        kwargs = pickle.load(f)

    # Load in some of the working state from sys.modules
    PAR = sys.modules["seisflows_parameters"]
    PATH = sys.modules["seisflows_paths"]
    system = sys.modules["seisflows_system"]

    # Configure the CPU-dependent logger which will log to stdout only
    # But mainsolver will log to the main log file as well
    if system.taskid == 0: 
        filename = PATH.LOGFILE
    else:
        filename = None
    config_logger(level=PAR.LOG_LEVEL, verbose=PAR.VERBOSE, filename=filename)

    # Get the actual function so we can evaluate it
    func = getattr(sys.modules[f"seisflows_{args.classname}"], args.funcname)

    # Evaluate the function with given keyword arguments
    func(**kwargs)


