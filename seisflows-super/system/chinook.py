#!/usr/bin/env python3
"""
Chinook is the University of Alaska Fairbanks (UAF) high performance computer,
operated by the Geophysical Institute's Research Computing Systems (RCS).
Chinook operates on a SLURM workload manager and therefore overloads the SLURM
System module. Chinook-specific parameters and functions are defined here.

Information on Chinook can be found here:
https://uaf-rcs.gitbook.io/uaf-rcs-hpc-docs/hpc
"""
import os
import sys
import logging

from seisflows3.config import custom_import, SeisFlowsPathsParameters

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Chinook(custom_import("system", "slurm")):
    """
    System interface for the University of Alaska HPC Chinook, which operates
    on a SLURM system.
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type partitions: dict
        :param partitions: Chinook has various partitions which each have their
            own number of cores per compute node, defined here
        """
        super().__init__()
        self.partitions = {"debug": 24, "t1small": 28, "t2small": 28,
                           "t1standard": 40, "t2standard": 40, "analysis": 28
                           }

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        sf.par("PARTITION", required=False, default="t1small", par_type=int,
               docstr="Name of partition on main cluster, available: "
                      "analysis, t1small, t2small, t1standard, t2standard, gpu")

        sf.par("MPIEXEC", required=False, default="srun", par_type=str,
               docstr="Function used to invoke parallel executables")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        if validate:
            self.required.validate()
        super().check(validate=False)

        assert(PAR.PARTITION in self.partitions.keys()), \
            f"Chinook partition must be in {self.partitions.keys()}"

        assert(PAR.NODESIZE == self.partitions[PAR.PARTITION]), \
            (f"PARTITION {PAR.PARTITION} is expected to have NODESIZE=" 
             f"{self.partitions[PAR.PARTITION]}, not current {PAR.NODESIZE}")

