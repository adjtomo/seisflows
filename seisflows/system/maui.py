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
from seisflows.system.slurm import Slurm
from seisflows.config import ROOT_DIR


class Maui(Slurm):
    """
    System interface for Maui, which operates on a SLURM system
    """
    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type partitions: dict
        :param partitions: Maui has various partitions which each have their
            own number of cores per compute node, defined here
        """
        super().__init__()

        self.required.par(
            "ACCOUNT", required=True, par_type=str,
            docstr="Maui account name to submit jobs under"
        )
        self.required.par(
            "NODESIZE", required=False, default=40, par_type=int,
            docstr="The number of cores per node defined by the Maui cluster. "
                   "Assumed to be 40 cores per node."
        )
        self.required.par(
            "MPIEXEC", required=False, default="srun", par_type=str,
            docstr="MPI call function used to invoke parallel executables, "
                   "defaults to 'srun'"
        )
        self.required.par(
            "CPUS_PER_TASK", required=False, default=1, par_type=int,
            docstr="Multiple CPUS per task allows for multithreading jobs"
        )
        self.required.par(
            "CLUSTER", required=False, default="maui", par_type=str,
            docstr="Name of main cluster for job submission. Available options: "
                   "'maui', 'maui_ancil', 'mahuika'. Note Mahuika untested"
        )
        self.required.par(
            "PARTITION", required=False, default="nesi_research",
            par_type=str, docstr="Name of cluster partition to submit job to"
        )
        self.required.par(
            "ANCIL_CLUSTER", required=False, default="maui_ancil", par_type=str,
            docstr="Ancillary cluster for pre- and post-processing tasks."
                   "Defaults to 'maui_ancil'")

        self.required.par(
            "ANCIL_PARTITION", required=False, default="nesi_prepost",
            par_type=str,
            docstr="Name of ancillary partition for prepost tasks. Defaults to"
                   "'nesi_prepost'"
        )
        self.required.par(
            "ANCIL_TASKTIME", required=False, default="null", par_type=float,
            docstr="Tasktime for prepost jobs submitted to ancillary nodes "
                   "matching 'ANCIL_CLUSTER' and 'ANCIL_PARTITION'"
        )

        self.partitions = {"nesi_research": 40}

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

        assert(self.par.NODESIZE == self.partitions[self.par.PARTITION]), \
            (f"PARTITION {self.par.PARTITION} is expected to have NODESIZE=" 
             f"{self.partitions[self.par.PARTITION]}, not current "
             f"{self.par.NODESIZE}")

        assert("SLURM_MEM_PER_CPU" in (self.par.ENVIRONS or "")), \
            ("Maui runs Slurm>=21 which enforces mutually exclusivity of Slurm "
             "memory environment variables SLURM_MEM_PER_CPU and "
             "SLURM_MEM_PER_NODE. Due to the cross-cluster nature of "
             "running SeisFlows3 on Maui, we must remove one env. variable. "
             "Please add 'SLURM_MEM_PER_CPU' to self.par.ENVIRONS.")

    def setup(self):
        """Inherits from workflow.system.workstation.Workstation"""
        self.setup()

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
                f"--account={self.par.ACCOUNT}",
                f"--cluster={self.par.ANCIL_CLUSTER}",
                f"--partition={self.par.ANCIL_PARTITION}",
                f"--job-name={self.par.TITLE}",
                f"--output={self.output_log}",
                f"--error={self.error_log}",
                f"--ntasks=1",
                f"--cpus-per-task=1",
                f"--time={self.par.WALLTIME:d}",
                f"{os.path.join(ROOT_DIR, 'scripts', 'submit')}",
                f"--output {self.path.OUTPUT}"
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
            _nodes = np.ceil(self.par.NPROC / float(self.par.NODESIZE))

            run_call = " ".join([
                "sbatch",
                f"{self.par.SLURMARGS or ''}",
                f"--account={self.par.ACCOUNT}",
                f"--job-name={self.par.TITLE}",
                f"--clusters={self.par.CLUSTER}",
                f"--partition={self.par.PARTITION}",
                f"--cpus-per-task={self.par.CPUS_PER_TASK}",
                f"--nodes={_nodes:d}",
                f"--ntasks={self.par.NPROC:d}",
                f"--time={self.par.TASKTIME:d}",
                f"--output={os.path.join(self.path.WORKDIR, 'logs', '%A_%a')}",
                f"--array=0-{self.par.NTASK-1 % self.par.NTASKMAX}",
                f"{os.path.join(ROOT_DIR, 'scripts', 'run')}",
                f"--output {self.path.OUTPUT}",
                f"--classname {classname}",
                f"--funcname {method}",
                f"--environment {self.par.ENVIRONS or ''}"
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
            f"{self.par.SLURMARGS or ''}",
            f"--account={self.par.ACCOUNT}",
            f"--job-name={self.par.TITLE}",
            f"--clusters={self.par.ANCIL_CLUSTER}",
            f"--partition={self.par.ANCIL_PARTITION}",
            f"--cpus-per-task={self.par.CPUS_PER_TASK}",
            f"--time={self.par.ANCIL_TASKTIME:d}",
            f"--output={os.path.join(self.path.WORKDIR, 'logs', '%A_%a')}",
            f"--array=0-{self.par.NTASK-1 % self.par.NTASKMAX}",
            f"{os.path.join(ROOT_DIR, 'scripts', 'run')}",
            f"--output {self.path.OUTPUT}",
            f"--classname {classname}",
            f"--funcname {method}",
            f"--environment {self.par.ENVIRONS or ''}"
        ])
        self.logger.debug(ancil_run_call)
        super().run(classname, method, single=False, run_call=ancil_run_call,
                    **kwargs)

    def taskid(self):
        """Inherits from seisflows.system.slurm.Slurm"""
        return self.taskid()

    def checkpoint(self, path, classname, method, kwargs):
        """Inherits from workflow.system.workstation.Workstation"""
        self.checkpoint(path=path, classname=classname, method=method,
                        kwargs=kwargs)



