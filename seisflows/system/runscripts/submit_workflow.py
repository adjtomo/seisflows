#!/usr/bin/env python3
"""
Only required when system==cluster (or any subclass of cluster)

This script is executes a MASTER job through job scheduler
(e.g., PBS, LSF, or SLURM) by running workflow.main() on the compute system.

.. note::
    Not to be called by the user. This is called when the user runs
    `seisflows submit` or `seisflows resume`. Internally this script is meant to
    be called by system.submit().

.. rubric::
    >> python submit --output ./OUTPUT
    OR
    >> sbatch submit --output ./OUTPUT
"""
import sys
import argparse

from seisflows.tools import unix
from seisflows.tools.config import load, config_logger


def parse_args():
    """
    Get command line arguments required for the submit script
    """
    parser = argparse.ArgumentParser("Run arguments for system submitted tasks")
    parser.add_argument("-o", "--output", type=str, nargs="?", required=True,
                        help="the SeisFlows output directory used to load the "
                             "active working state from inside the compute node"
                        )

    return parser.parse_args()


if __name__ == '__main__':
    """
    Submit workflow.main() as a MASTER JOB on the compute system
    """
    args = parse_args()

    # Load the currently active working state
    unix.cd(args.output)
    load(args.output)

    # Ensure that the two main modules are loaded
    workflow = sys.modules["seisflows_workflow"]
    system = sys.modules["seisflows_system"]

    # Set up logging on the compute system to print to stdout only
    PAR = sys.modules["seisflows_parameters"]
    PATH = sys.modules["seisflows_paths"]
    config_logger(level=PAR.LOG_LEVEL, verbose=PAR.VERBOSE) 

    # Execute MASTER JOB as workflow.main()
    workflow.main()

