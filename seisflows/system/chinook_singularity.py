#!/usr/bin/env python3
"""
For running SeisFlows on Chinook through Singularity. This class is a
kludge to get research going as Chinook is scheduled for OS upgrade in the
coming month(s) which will cause this entire kludge to become obsolete. This
code is therefore unpolished and only meant to get things working in a
band-aid/barebones manner.
"""
import os
from datetime import timedelta
from seisflows.system.singularity import Singularity


class ChinookSingularity(Singularity):
    """
    System Chinook Singularity
    --------------------------
    University of Alaska Fairbanks HPC Chinook, SLURM based system running
    the Singularity container software

    Parameters
    ----------
    :type partition: str
    :param partition: Chinook has various partitions which each have their
        own number of cores per compute node. Available are: analysis, t1small,
        t2small, t1standard, t2standard, gpu

    Paths
    -----

    ***
    """
    __doc__ = Singularity.__doc__ + __doc__


    def __init__(self, partition="t1small", slurm_args="", **kwargs):
        """Chinook init"""
        super().__init__(**kwargs)

        self.partition = partition
        self.slurm_args = slurm_args  # currently NOT used

        self._partitions = {"debug": 24, "t1small": 28, "t2small": 28,
                           "t1standard": 40, "t2standard": 40, "analysis": 28
                           }
        self._tasktime = str(timedelta(minutes=self.tasktime))


    @property
    def node_size(self):
        """Defines the node size of a given cluster partition. This is a hard
        set number defined by the system architecture"""
        return self._partitions[self.partition]

    @property
    def run_call_header(self):
        """
        The run call defines the SBATCH header which is used to run tasks during
        an executing workflow. Like the submit call its arguments are dictated
        by the given system. Run calls are modified and called by the `run`
        function

        :rtype: str
        :return: the system-dependent portion of a run call
        """
        _output = os.path.relpath(os.path.join(self.path.log_files, '%A_%a'))
        _call = "\n".join([
            f"#SBATCH --job-name={self.title}",
            f"#SBATCH --ntasks={self.nproc:d}",
            f"#SBATCH --tasks-per-node={self.node_size}",
            f"#SBATCH --time={self._tasktime}",
            f"#SBATCH --output={_output}",
            f"#SBATCH --array=0-{self.ntask-1 % self.ntask_max}",
            f"{super().run_call_header}",  # Singularity wrapper
        ])

        return _call
