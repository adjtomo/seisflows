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
from seisflows.system.slurm import Slurm


class Chinook(Slurm):
    """
    System Chinook
    --------------
    University of Alaska Fairbanks HPC Chinook, SLURM based system

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
    __doc__ = Slurm.__doc__ + __doc__


    def __init__(self, partition="t1small", **kwargs):
        """Chinook init"""
        super().__init__(**kwargs)

        self.partition = partition

        self._partitions = {"debug": 24, "t1small": 28, "t2small": 28,
                            "t1standard": 40, "t2standard": 40, "analysis": 28
                            }

    @property
    def submit_call_header(self):
        """
        The submit call defines the SBATCH header which is used to submit a
        workflow task list to the system. It is usually dictated by the
        system's required parameters, such as account names and partitions.
        Submit calls are modified and called by the `submit` function.

        .. note::
            Force the partition to run on 'Debug' since we are only running a 
            single core job. This may run into walltime problems but for now
            it seems more economincal to enforce this

        :rtype: str
        :return: the system-dependent portion of a submit call
        """
        _call = " ".join([
            f"sbatch",
            f"--job-name={self.title}",
            f"--output={self.path.output_log}",
            f"--error={self.path.output_log}",
            f"--ntasks=1",
            f"--partition=debug",
            # f"--partition={self.partition}",
            f"--time={self._walltime}"
        ])
        return _call

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
        _call = " ".join([
            f"sbatch",
            f"{self.slurm_args or ''}",
            f"--job-name={self.title}",
            f"--ntasks={self.nproc:d}",
            f"--partition={self.partition}",
            f"--tasks-per-node={self.node_size}",
            f"--time={self._tasktime}",
            f"--output={os.path.join(self.path.log_files, '%A_%a')}",
            f"--array=0-{self.ntask-1 % self.ntask_max}",
            f"--parsable"
        ])
        return _call
