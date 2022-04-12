#!/usr/bin/env python3
"""
Only required when system==cluster (or any subclass of cluster)

This script is a wrapper for running tasks on systems during an active workflow.
Acts as a Python script to submit certain SeisFlows3 functions or tasks to a
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

from seisflows3.config import load, config_logger


def parse_args():
    """
    Get command line arguments
    """
    parser = argparse.ArgumentParser("Run arguments for system submitted tasks")
    parser.add_argument("-o", "--output", type=str, nargs="?", required=True,
                        help="the SeisFlows3 output directory used to load the "
                             "active working state from inside the compute node"
                        )
    parser.add_argument("-c", "--classname", type=str, nargs="?", required=True,
                        help="the SeisFlows3 class from within which the "
                             "desired function is defined. Available options "
                             "are defined in seisflows3.config.NAMES"
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
    Exports comma delimited list of environment variables

    e.g. VARNAME1=value1,VARNAME2=value2

    :type myenv: str
    :param myevn: the system environment to take variables from
    """
    for item in myenv.split(","):
        os.environ.update([item.split("=")])


if __name__ == '__main__':
    """ 
    Runs task within a currently executing workflow 
    """
    args = parse_args()

    if args.environment:
        export(args.environment)

    # Load the last checkpointed working state from the 'seisflows_?.p` files
    load(args.output)

    # Load keyword arguments required by this function
    # Files will be something like: 'solver_eval_func.p'
    kwargs_fid = f"{args.classname}_{args.funcname}.p"
    kwargs_path = os.path.join(args.output, "kwargs", kwargs_fid)
    with open(kwargs_path, "rb") as f:
        kwargs = pickle.load(f)

    # Configure the CPU-dependent logger
    PAR = sys.modules["seisflows_parameters"]
    config_logger(level=PAR.LOG_LEVEL, verbose=PAR.VERBOSE)

    # Get the actual function so we can evaluate it
    func = getattr(sys.modules[f"seisflows_{args.classname}"], args.funcname)

    # Evaluate the function with given keyword arguments
    func(**kwargs)


