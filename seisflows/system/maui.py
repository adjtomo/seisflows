#!/usr/bin/env python3
"""
Maui is a New Zealand eScience Infrastructure (NeSI) high performance computer.
Maui operates on a SLURM workload manager and therefore overloads the SLURM
System module. Maui-specific parameters and functions are defined here.

Information on Maui can be found here:
https://support.nesi.org.nz/hc/en-gb/articles/360000163695-M%C4%81ui

.. note::
    Python and conda capabilities are NOT accessible from Maui, these
    capabilities have been shifted onto a separate cluster: Maui ancil
    This subclass therefore moves all Python dependent capabilities
    (i.e., SeisFlows3, Pyatoa) onto the ancilary cluster.

    See also: https://support.nesi.org.nz/hc/en-gb/articles/\
                                          360000203776-M%C4%81ui-Ancillary-Nodes

"""
import os
from seisflows.system.slurm import Slurm


class Maui(Slurm):
    """
    System Maui
    -----------
    New Zealand Maui-specfic modifications to base SLURM system

    Parameters
    ----------
    :type account: str
    :param account: Maui account to submit jobs under, will be used for the
        '--account' sbatch argument
    :type cpus_per_task: int
    :param cpus_per_task: allow for multiple cpus per task, i.e,.
        multithreaded jobs
    :type cluster: str
    :param cluster: cluster to submit jobs to. Available are Maui and
        Mahuika
    :type partition: str
    :param partition: partition of the cluster to submit jobs to.
    :type ancil_cluster: str
    :param ancil_cluster: name of the ancilary cluster used for pre-
        post-processing tasks.
    :type ancil_partition: name of the partition of the ancilary cluster
    :type ancil_tasktime: int
    :param ancil_tasktime: Tasktime in minutes for pre and post-processing
        jobs submitted to Maui ancil.

    Paths
    -----
    ***
    """
    __doc__ = Slurm.__doc__ + __doc__

    def __init__(self, account=None, cpus_per_task=1, cluster="maui",
                 partition="nesi_research", ancil_cluster="maui_ancil",
                 ancil_partition="nesi_prepost", ancil_tasktime=1, **kwargs):
        """Maui init"""
        super().__init__(**kwargs)

        self.account = account
        self.cluster = cluster
        self.partition = partition
        self.cpus_per_task = cpus_per_task
        self.ancil_cluster = ancil_cluster
        self.ancil_partition = ancil_partition
        self.ancil_tasktime = ancil_tasktime
        if self.environs and "SLURM_MEM_PER_CPU" not in self.environs:
            self.environs = f"{self.environs},SLURM_MEM_PER_CPU"
        else:
            self.environs = "SLURM_MEM_PER_CPU"

        self._partitions = {"nesi_research": 40}
        self._available_clusters = ["maui", "mahuika"]

    def check(self):
        """
        Checks parameters and paths
        """
        super().check()

        assert("SLURM_MEM_PER_CPU" in (self.environs or "")), \
            ("Maui runs Slurm>=21 which enforces mutually exclusivity of Slurm "
             "memory environment variables SLURM_MEM_PER_CPU and "
             "SLURM_MEM_PER_NODE. Due to the cross-cluster nature of "
             "running SeisFlows on Maui, we must remove one env. variable. "
             "Please add 'SLURM_MEM_PER_CPU' to self.par.ENVIRONS.")

        assert(self.cluster in self._available_clusters), (
            f"System 'Maui' parameter cluster must be in "
            f"{self._available_clusters}"
            )

        assert(self.account), f"System 'Maui' requires parameter 'account'"

    @property
    def submit_call_header(self):
        """
        The submit call defines the SBATCH header which is used to submit a
        workflow task list to the system. It is usually dictated by the
        system's required parameters, such as account names and partitions.
        Submit calls are modified and called by the `submit` function.

        .. note::
            The master job must be run on `maui_ancil` because Maui does
            not have the ability to run the command "sacct", nor can it
            not have the ability to run the command "sacct", nor can it
            use the Conda environment that has been set by Ancil

        .. note::
            We do not place SLURMARGS into the sbatch command to avoid the
            export=None which will not propagate the conda environment

        :rtype: str
        :return: the system-dependent portion of a submit call
        """
        _call = " ".join([
            f"sbatch",
            f"--account={self.account}",
            f"--cluster={self.ancil_cluster}",
            f"--partition={self.ancil_partition}",
            f"--job-name={self.title}",
            f"--output={self.path.output_log}",
            f"--error={self.path.output_log}",
            f"--ntasks=1",
            f"--cpus-per-task=1",
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
             f"--account={self.account}",
             f"--job-name={self.title}",
             f"--clusters={self.cluster}",
             f"--partition={self.partition}",
             f"--cpus-per-task={self.cpus_per_task}",
             f"--nodes={self.nodes:d}",
             f"--ntasks={self.nproc:d}",
             f"--time={self._tasktime}",
             f"--output={os.path.join(self.path.log_files, '%A_%a')}",
             f"--array=0-{self.ntask-1 % self.ntask_max}",
             f"--parsable"
        ])
        return _call

    @property
    def ancil_run_call_header(self):
        """
        A modified form of `run_call` which is used to run jobs on the Ancil
        pre/postprocessing cluster of Maui. This is used to run Pyaflowa jobs
        which require the Conda environment active on Maui Ancil.
        """
        _call = " ".join([
             f"sbatch",
             f"{self.slurm_args or ''}",
             f"--account={self.account}",
             f"--job-name={self.title}",
             f"--clusters={self.ancil_cluster}",
             f"--partition={self.ancil_partition}",
             f"--cpus-per-task={self.cpus_per_task}",
             f"--nodes={self.nodes:d}",
             f"--ntasks={self.nproc:d}",
             f"--time={self.ancil_tasktime:d}",
             f"--output={os.path.join(self.path.log_files, '%A_%a')}",
             f"--array=0-{self.ntask-1 % self.ntask_max}"
        ])
        return _call
