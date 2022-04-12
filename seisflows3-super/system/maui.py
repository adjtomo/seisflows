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
import sys
import math
import logging
from seisflows3.config import custom_import, SeisFlowsPathsParameters, ROOT_DIR

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Maui(custom_import("system", "slurm")):
    """
    System interface for Maui, which operates on a SLURM system
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type partitions: dict
        :param partitions: Maui has various partitions which each have their
            own number of cores per compute node, defined here
        """
        super().__init__()
        self.partitions = {"nesi_research": 40}

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        sf.par("ACCOUNT", required=True, par_type=str,
               docstr="The account name to submit jobs under")

        sf.par("CPUS_PER_TASK", required=False, default=1, par_type=int,
               docstr="Multiple CPUS per task allows for multithreading jobs")

        sf.par("CLUSTER", required=False, default="maui", par_type=str,
               docstr="Name of main cluster for parallel job submission")

        sf.par("PARTITION", required=False, default="nesi_research",
               par_type=str, docstr="Name of partition on main cluster")

        sf.par("ANCIL_CLUSTER", required=False, default="maui_ancil",
               par_type=str,
               docstr="Name of ancillary cluster for prepost tasks")

        sf.par("ANCIL_PARTITION", required=False, default="nesi_prepost",
               par_type=str,
               docstr="Name of ancillary partition for prepost tasks")

        sf.par("ANCIL_TASKTIME", required=False, default="null", par_type=float,
               docstr="Tasktime for prepost jobs on ancillary nodes")

        sf.par("NODESIZE", required=False, default=40, par_type=int,
               docstr="The number of cores per node defined by the system")

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

        assert(PAR.NODESIZE == self.partitions[PAR.PARTITION]), \
            (f"PARTITION {PAR.PARTITION} is expected to have NODESIZE=" 
             f"{self.partitions[PAR.PARTITION]}, not current {PAR.NODESIZE}")

    def submit(self):
        """
        Submits master job workflow to maui_ancil cluster as a single-core
        process

        .. note::
            The master job must be run on maui_ancil because Maui does
            not have the ability to run the command "sacct"
        """
        maui_submit_call = " ".join([
            f"sbatch {PAR.SLURMARGS or ''}",
            f"--account={PAR.ACCOUNT}",
            f"--cluster={PAR.ANCIL_CLUSTER}",
            f"--partition={PAR.ANCIL_PARTITION}",
            f"--job-name={PAR.TITLE}",
            f"--output=output-%A.txt",
            f"--error=error-%A.txt",
            f"--ntasks=1",
            f"--cpus-per-task=1",
            f"--time={PAR.WALLTIME:d}",
            f"{os.path.join(ROOT_DIR, 'scripts', 'submit')}",
            f"--output {PATH.OUTPUT}"
        ])
        self.logger.debug(maui_submit_call)
        super().submit(maui_submit_call)

    def run(self, classname, method, single=False, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on a SLURM
        cluster. Executes classname.method(*args, **kwargs) `NTASK` times,
        each time on `NPROC` CPU cores

        :type classname: str
        :param classname: the class to run
        :type method: str
        :param method: the method from the given `classname` to run
        :type scale_tasktime: int
        :param scale_tasktime: a way to get over the hard-set tasktime, because
            some tasks take longer (e.g. smoothing), but you don't want these
            to set the tasktimes for all other tasks. This lets you scale the
            time of specific tasks by PAR.TASKTIME * scale_tasktime
        """
        maui_run_call = " ".join([
            "sbatch",
            f"{PAR.SLURMARGS or ''}",
            f"--account={PAR.ACCOUNT}",
            f"--job-name={PAR.TITLE}",
            f"--clusters={PAR.CLUSTER}",
            f"--partition={PAR.PARTITION}",
            f"--cpus-per-task={PAR.CPUS_PER_TASK}",
            f"--nodes={math.ceil(PAR.NPROC / float(PAR.NODESIZE)):d}",
            f"--ntasks={PAR.NPROC:d}",
            f"--time={PAR.TASKTIME:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'logs', '%A_%a')}",
            f"--array=0-{PAR.NTASK-1 % PAR.NTASKMAX}",
            f"{os.path.join(ROOT_DIR, 'scripts', 'run')}",
            f"--output {PATH.OUTPUT}",
            f"--classname {classname}",
            f"--funcname {method}",
            f"--environment {PAR.ENVIRONS}"
        ])
        self.logger.debug(maui_run_call)
        super().run(classname, method, single, run_call=maui_run_call, **kwargs)

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
            f"{PAR.SLURMARGS or ''}",
            f"--account={PAR.ACCOUNT}",
            f"--job-name={PAR.TITLE}",
            f"--clusters={PAR.ANCIL_CLUSTER}",
            f"--partition={PAR.ANCIL_PARTITION}",
            f"--cpus-per-task={PAR.CPUS_PER_TASK}",
            f"--time={PAR.ANCIL_TASKTIME:d}",
            f"--output={os.path.join(PATH.WORKDIR, 'logs', '%A_%a')}",
            f"--array=0-{PAR.NTASK-1 % PAR.NTASKMAX}",
            f"{os.path.join(ROOT_DIR, 'scripts', 'run')}",
            f"--output {PATH.OUTPUT}",
            f"--classname {classname}",
            f"--funcname {method}",
            f"--environment {PAR.ENVIRONS}"
        ])
        self.logger.debug(ancil_run_call)
        super().run(classname, method, single=False, run_call=ancil_run_call,
                    **kwargs)



    

