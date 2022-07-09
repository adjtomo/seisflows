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
import numpy as np
from seisflows import logger
from seisflows.system.slurm import Slurm
from seisflows.config import ROOT_DIR


class Maui(Slurm):
    """
    System interface for Maui, which operates on a SLURM system
    """
    def __init__(self, account=None, cpus_per_task=1, cluster="maui",
                 partition="nesi_research", ancil_cluster="maui_ancil",
                 ancil_partition="nesi_prepost", ancil_tasktime=1, **kwargs):
        """
        Maui parameters

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
        """
        super().__init__(**kwargs)

        self.account = account
        self.cluster = cluster
        self.partition = partition
        self.cpus_per_task = cpus_per_task
        self.ancil_cluster = ancil_cluster
        self.ancil_partition = ancil_partition
        self.ancil_tasktime = ancil_tasktime

        self._partitions = {"nesi_research": 40}
        self.node_size = self._partitions[self.partition]

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

        assert("SLURM_MEM_PER_CPU" in (self.environs or "")), \
            ("Maui runs Slurm>=21 which enforces mutually exclusivity of Slurm "
             "memory environment variables SLURM_MEM_PER_CPU and "
             "SLURM_MEM_PER_NODE. Due to the cross-cluster nature of "
             "running SeisFlows3 on Maui, we must remove one env. variable. "
             "Please add 'SLURM_MEM_PER_CPU' to self.par.ENVIRONS.")

    def submit(self, submit_call=None):
        """
        Submits master job workflow to maui_ancil cluster as a single-core
        process

        .. note::
            The master job must be run on maui_ancil because Maui does
            not have the ability to run the command "sacct", nor can it
            use the Conda environment that has been set by Ancil

        .. note::
            We do not place SLURMARGS into the sbatch command to avoid the
            export=None which will not propagate the conda environment

        :type submit_call: str
        :param submit_call: SBATCH command line call to submit workflow.main()
            to the system. If None, will generate one on the fly with
            user-defined parameters
        """
        if submit_call is None:
            submit_call = " ".join([
                f"sbatch",
                f"--account={self.account}",
                f"--cluster={self.ancil_cluster}",
                f"--partition={self.ancil_partition}",
                f"--job-name={self.title}",
                f"--output={self.path_output_log}",
                f"--error={self.path_error_log}",
                f"--ntasks=1",
                f"--cpus-per-task=1",
                f"--time={self.walltime:d}",
                f"{os.path.join(ROOT_DIR, 'system', 'runscripts', 'submit')}",
                f"--output {self.path_output}"
            ])

        super().submit(submit_call=submit_call)

    def run(self, classname, method, single=False, run_call=None, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on a SLURM
        cluster. Executes classname.method(*args, **kwargs) `NTASK` times,
        each time on `NPROC` CPU cores

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        :type single: bool
        :param single: run a single-process, non-parallel task, such as
            smoothing the gradient, which only needs to be run by once.
            This will change how the job array and the number of tasks is
            defined, such that the job is submitted as a single-core job to
            the system.
        :type run_call: str
        :param run_call: SBATCH command line run call to be submitted to the
            system. If None, will generate one on the fly with user-defined
            parameters
        """
        if run_call is None:
            # Calculate requested number of nodes based on requested proc count
            _nodes = np.ceil(self.nproc / float(self.node_size))
            _nodes = _nodes.astype(int)

            run_call = " ".join([
                "sbatch",
                f"{self.slurm_args or ''}",
                f"--account={self.account}",
                f"--job-name={self.title}",
                f"--clusters={self.cluster}",
                f"--partition={self.partition}",
                f"--cpus-per-task={self.cpus_per_task}",
                f"--nodes={_nodes:d}",
                f"--ntasks={self.nproc:d}",
                f"--time={self.tasktime:d}",
                f"--output={os.path.join(self.path_log_files, '%A_%a')}",
                f"--array=0-{self.ntask-1 % self.ntask_max}",
                f"{os.path.join(ROOT_DIR, 'system', 'runscripts', 'run')}",
                f"--output {self.path_output}",
                f"--classname {classname}",
                f"--funcname {method}",
                f"--environment {self.environs or ''}"
            ])

        super().run(classname, method, single, run_call=run_call, **kwargs)

    def run_ancil(self, classname, method, **kwargs):
        """
        Runs prepost jobs on Maui ancil, the ancilary cluster which contains
        the conda and Python capabilities for Maui.

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        """
        ancil_run_call = " ".join([
            "sbatch",
            f"{self.slurm_args or ''}",
            f"--account={self.account}",
            f"--job-name={self.title}",
            f"--clusters={self.ancil_cluster}",
            f"--partition={self.ancil_partition}",
            f"--cpus-per-task={self.cpus_per_task}",
            f"--time={self.ancil_tasktime:d}",
            f"--output={os.path.join(self.path_log_files, '%A_%a')}",
            f"--array=0-{self.ntask-1 % self.ntask_max}",
            f"{os.path.join(ROOT_DIR, 'system', 'runscripts', 'run')}",
            f"--output {self.path_output}",
            f"--classname {classname}",
            f"--funcname {method}",
            f"--environment {self.environs or ''}"
        ])
        logger.debug(ancil_run_call)
        super().run(classname, method, single=False, run_call=ancil_run_call,
                    **kwargs)
