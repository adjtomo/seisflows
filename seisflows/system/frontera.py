#!/usr/bin/env python3
"""
Frontera is one of the Texas Advanced Computing Center (TACC) HPCs.
https://frontera-portal.tacc.utexas.edu/
"""
import os
from seisflows.system.slurm import Slurm


class Frontera(Slurm):
    """
    System Frontera
    --------------
    Texas Advanced Computing Center HPC Frontera, SLURM based system

    Parameters
    ----------
    :type partition: str
    :param partition: Chinook has various partitions which each have their
        own number of cores per compute node. Available are: small, normal,
        large, development, flex
    :type allocation: str
    :param allocation: Name of allocation/project on the Frontera system.
        Required if you have more than one active allocation.

    Paths
    -----

    ***
    """
    def __init__(self, partition="small", allocation=None, **kwargs):
        """Frontera init"""
        super().__init__(**kwargs)

        self.partition = partition
        self.allocation = allocation
        self.mpiexec = "ibrun"

        # TODO find out the cores-per-node values for these partitions
        self.partitions = {"small": None, "normal": None, "large": None,
                           "development": None, "flex": None}


    @property
    def submit_call_header(self):
        """
        The submit call defines the SBATCH header which is used to submit a
        workflow task list to the system. It is usually dictated by the
        system's required parameters, such as account names and partitions.
        Submit calls are modified and called by the `submit` function.

        :rtype: str
        :return: the system-dependent portion of a submit call
        """
        _call = " ".join([
            f"sbatch",
            f"--job-name={self.title}",  # -J
            f"--partition={self.partition}",  # -p
            f"--output={self.path.output_log}",  # -o
            f"--error={self.path.output_log}",
            f"--nodes=1",  # -N
            f"--ntasks=1",  # -n
            f"--time={self._walltime}"  # -t
        ])
        if self.allocation is not None:
            _call = f"{_call} --allocation={self.allocation}"
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
            f"--partition={self.partition}",
            f"--output={os.path.join(self.path.log_files, '%A_%a')}",
            f"--ntasks={self.nproc:d}",
            f"--nodes={self.nodes}",
            f"--array=0-{self.ntask - 1 % self.ntask_max}",
            f"--time={self._tasktime}",
            f"--parsable"
        ])
        if self.allocation is not None:
            _call = f"{_call} --allocation={self.allocation}"
        return _call
