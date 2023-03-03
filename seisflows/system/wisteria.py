#!/usr/bin/env python3
"""
Wisteria is the University of Tokyo Fujitsu brand high performance computer.
Wisteria runs on the Fujitsu/PJM job scheduler.

.. notes::

    - Wisteria has two node gruops, Odyssey (compute nodes) and Aquarius 
      (data/learning nodes w/ GPU)
    - Odyssey has 7680 nodes with 48 cores/node
    - Aquarius has 45 nodes with 36 cores/node
"""
import os
from seisflows.system.fujitsu import Fujitsu


class Wisteria(Fujitsu):
    """
    System Wisteria
    ---------------
    University of Tokyo HPC Wisteria, running Fujitsu job scheduler

    Parameters
    ----------

    Paths
    -----

    ***
    """
    __doc__ = Slurm.__doc__ + __doc__


    def __init__(self, rscgrp="t1small", **kwargs):
        """Chinook init"""
        super().__init__(**kwargs)

        self.mpiexec = mpiexec
        self.partition = partition
        self.submit_to = submit_to or self.partition

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
            f"--partition={self.submit_to}",
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
            # f"--tasks-per-node={self.node_size}",  # actually not required?
            f"--time={self._tasktime}",
            f"--output={os.path.join(self.path.log_files, '%A_%a')}",
            f"--array=0-{self.ntask-1 % self.ntask_max}",
            f"--parsable"
        ])
        return _call
