#!/usr/bin/env python3
"""
Only required when system==cluster (or any subclass of cluster)
This script is used to execute a MASTER job on system. It is essentially the
same as `system.workstation.Workstation.submit()`, except it can be called as
a script using subprocess, allowing it to be submitted to e.g., a job scheduler

.. note::
    Not to be called by the user. This is called when the user runs
    `seisflows submit` or `seisflows resume`. Internally this script is meant to
    be called by system.submit().

.. rubric::
    >> python submit -w ./ -p parameters.yaml
    OR
    >> sbatch submit -w ./ -p parameters.yaml
"""
import os
import argparse
from seisflows.tools.config import import_seisflows


def parse_args():
    """
    Get command line arguments required for the submit script
    """
    parser = argparse.ArgumentParser("Run arguments for system submitted tasks")
    parser.add_argument("-w", "--workdir", type=str, nargs="?", required=True,
                        default=os.getcwd(), help="SeisFlows working directory")
    parser.add_argument("-p", "--parameter_file", type=str, nargs="?",
                        required=True, default="parameters.yaml",
                        help="SeisFlows parameter file")

    return parser.parse_args()


if __name__ == '__main__':
    """
    Submit workflow.main() as a MASTER JOB on the compute system
    """
    args = parse_args()
    workflow = import_seisflows(workdir=args.workdir,
                                parameter_file=args.parameter_file,
                                stream_handler=False)
    workflow.check()
    workflow.setup()
    workflow.run()
