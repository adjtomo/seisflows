#!/usr/bin/env python3
"""
Frontera is one of the Texas Advanced Computing Center (TACC) HPCs.
https://frontera-portal.tacc.utexas.edu/

TODO we may need to include or create a "singularity" class or run script which
runs jobs through singularity
"""
import os
import numpy as np
from seisflows.config import ROOT_DIR
from seisflows.system.slurm import Slurm


class Frontera(Slurm):
    """
    System interface for TACC Frontera based on SLURM workload manager
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
            "PARTITION", required=False, default="small", par_type=str,
            docstr="Name of partition on main cluster"
        )
        self.required.par(
            "ALLOCATION", required=False, default="", par_type=str,
            docstr="Name of allocation/project on the Frontera system. "
                   "Required if you have more than one active allocation."
        )
        self.required.par(
            "MPIEXEC", required=False, default="ibrun", par_type=str,
            docstr="Function used to invoke parallel executables. Defaults to"
                   "'ibrun' based on TACC user manual.")

        # TODO find out the cores-per-node values for these partitions
        # self.partitions = {"small":, "normal":, "large":, "development:"
        #                    "flex":}

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

    def submit(self, submit_call=None):
        """
        Submits workflow as a serial job on the TACC partition 'small'.

        .. note::
            The SBATCH commands can either be short or full length. TACC's
            start up guide uses short length keys so that's what we do here, but
            their long names can be substituted

        :type submit_call: str
        :param submit_call: SBATCH command line call to submit workflow.main()
            to the system. If None, will generate one on the fly with
            user-defined parameters
        """
        if submit_call is None:
            submit_call = " ".join([
                "sbatch",
                f"{self.par.SLURMARGS or ''}",
                f"-J {self.par.TITLE}",  # job name
                f"-O {self.output_log}",  # stdout output file
                f"-E {self.error_log}",  # stderr error file
                f"-P {self.par.PARTITION}",  # queue/partition name
                f"-A {self.par.ALLOCATION}",  # project/allocation name
                f"-N 1",  # total number of nodes requested
                f"-n 1",  # number of mpi tasks
                f"-t {self.par.WALLTIME}",  # job walltime
                f"{os.path.join(ROOT_DIR, 'scripts', 'submit')}",
                f"--output {self.path.OUTPUT}"
            ])
        super().submit(submit_call=submit_call)

    def run(self, classname, method, single=False, run_call=None, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on a SLURM
        cluster.

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
                f"-J {self.par.TITLE}",  # job name
                f"-O {self.output_log}",  # stdout output file
                f"-E {self.error_log}",  # stderr error file
                f"-P {self.par.PARTITION}",  # queue/partition name
                f"-A {self.par.ALLOCATION}",  # project/allocation name
                f"-N {_nodes}",  # total number of nodes requested
                f"-n {self.par.NPROC}",  # number of mpi tasks
                f"-t {self.par.WALLTIME}",  # job walltime
                f"{os.path.join(ROOT_DIR, 'scripts', 'run')}",
                f"--output {self.path.OUTPUT}"
                f"--classname {classname}",
                f"--funcname {method}",
                f"--environment {self.par.ENVIRONS or ''}"
            ])

        super().run(classname, method, single, run_call=run_call, **kwargs)

    def taskid(self):
        """Inherits from seisflows.system.slurm.Slurm"""
        return self.taskid()

    def checkpoint(self, path, classname, method, kwargs):
        """Inherits from workflow.system.workstation.Workstation"""
        self.checkpoint(path=path, classname=classname, method=method,
                        kwargs=kwargs)
