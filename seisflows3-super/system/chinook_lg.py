#!/usr/bin/env python
"""
This is the subclass seisflows.system.chinook_lg
This class provides the core utilities interaction with HPC systems which must
be overloaded by subclasses
"""
import os
import sys
import math
import time
from glob import glob
from subprocess import check_output, call, CalledProcessError

from seisflows3.tools import unix
from seisflows3.tools.wrappers import call, findpath
from seisflows3.config import custom_import

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class ChinookLg(custom_import('system', 'slurm_lg')):
    """
    System interface for the University of Alaska cluster, Chinook.

    Pre/postprocessing tasks must be run on separate partitions
    """
    def check(self):
        """
        Checks parameters and paths
        """
        # Run SlurmLG checks first
        super().check()

        # RCS Nodesize is set at 24 or 28 cores per node
        if PAR.NODESIZE != 24:
            print("Chinook must have a nodesize of 24, overwriting user set")
            setattr(PAR, "NODESIZE", 24)

        # How to invoke executables
        if "MPIEXEC" not in PAR:
            setattr(PAR, "MPIEXEC", "mpiexec")

        # Specific partition of the main cluster
        if "MAIN_PARTITION" not in PAR:
            setattr(PAR, "MAIN_PARTITION", "t1small")

        # Specific partition of ancilary cluster
        if "ANCIL_PARTITION" not in PAR:
            setattr(PAR, "ANCIL_PARTITION", "analysis")

        # If preprocessing tasks are much less than PAR.TASKTIME, it can be
        # useful to manually set a shorter tasktime
        if "ANCIL_TASKTIME" not in PAR:
            setattr(PAR, "ANCIL_TASKTIME", PAR.TASKTIME)

        # If number of nodes not given, automatically calculate.
        # if the "nomultithread" hint is given, the number of nodes will need 
        # to be manually set
        if "NODES" not in PAR:
            setattr(PAR, "NODES", math.ceil(PAR.NPROC/float(PAR.NODESIZE)))

        # Optional additional SLURM arguments
        if "SLURMARGS" not in PAR:
            setattr(PAR, "SLURMARGS", "")

    def submit(self, workflow):
        """
        Overwrites seisflows.workflow.slurm_lg.submit()

        Submits master job workflow to main cluster
        """
        output_log, error_log = self.setup()
        workflow.checkpoint()

        # Submit to maui_ancil
        submit_call = " ".join([
            f"sbatch {PAR.SLURMARGS}",
            f"--partition={PAR.MAIN_PARTITION}",
            f"--job-name=main_{PAR.TITLE}",  # main_ prefix means master
            f"--output={output_log}-%A.log",
            f"--error={error_log}-%A.log",
            f"--nodes=1",
            f"--ntasks-per-node={PAR.NODESIZE}",
            f"--time={PAR.WALLTIME:d}",
            os.path.join(findpath("seisflows.system"), "wrappers", "submit"),
            PATH.OUTPUT
        ])
        call(submit_call)
