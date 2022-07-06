#!/usr/bin/env python3
"""
Chinook is the University of Alaska Fairbanks (UAF) high performance computer,
operated by the Geophysical Institute's Research Computing Systems (RCS).
Chinook operates on a SLURM workload manager and therefore overloads the SLURM
System module. Chinook-specific parameters and functions are defined here.

Information on Chinook can be found here:
https://uaf-rcs.gitbook.io/uaf-rcs-hpc-docs/hpc
"""
from seisflows.system.slurm import Slurm


class Chinook(Slurm):
    """
    System interface for the University of Alaska HPC Chinook, which operates
    on a SLURM system.
    """
    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type partitions: dict
        :param partitions: Chinook has various partitions which each have their
            own number of cores per compute node, defined here
        """
        super().__init__()

        self.required.par(
            "PARTITION", required=False, default="t1small", par_type=int,
            docstr="Name of partition on main cluster, available: "
                   "analysis, t1small, t2small, t1standard, t2standard, gpu")

        self.required.par(
            "MPIEXEC", required=False, default="srun", par_type=str,
            docstr="Function used to invoke parallel executables")

        self.partitions = {"debug": 24, "t1small": 28, "t2small": 28,
                           "t1standard": 40, "t2standard": 40, "analysis": 28
                           }

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

        assert(self.par.PARTITION in self.partitions.keys()), \
            f"Chinook partition must be in {self.partitions.keys()}"

        assert(self.par.NODESIZE == self.partitions[self.par.PARTITION]), \
            (f"PARTITION {self.par.PARTITION} is expected to have NODESIZE=" 
             f"{self.partitions[self.par.PARTITION]}, not current "
             f"{self.par.NODESIZE}")

