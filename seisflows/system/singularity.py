#!/usr/bin/env python3
"""
A Cluster-adjacent base class that provides core utilities for interactions
with HPC systems running Singularity. Must be overloaded by subclasses defined
for specific workload managers / clusters.

The `Singularity` class was written for clusters running SeisFlows through
Docker containers using Singularity. The reason for writing a separate class
is because Docker containers do not have access to the workload manager (e.,g
SLURM/sbatch) and therefore we cannot run job submission calls directly from
the Python environment. Instead, each time a job must be submitted to the
Cluster, the User must manually submit.

.. note::
    To users looking to run SeisFlows directly via their Cluster Conda
    environment, look at the `Cluster` class and its workload manager-specific
    sub-classes
"""
import os
import subprocess
from seisflows import logger, ROOT_DIR
from seisflows.tools import msg
from seisflows.tools.unix import nproc, cp
from seisflows.tools.config import pickle_function_list, import_seisflows
from seisflows.system.workstation import Workstation
from seisflows.system.slurm import modify_run_call_single_proc


class Singularity(Workstation):
    """
    Singularity System
    ------------------
    HPC interfacing through Docker/Singularity containers

    Parameters
    ----------
    :type title: str
    :param title: The name used to submit jobs to the system, defaults
        to the name of the current working directory
    :type mpiexec: str
    :param mpiexec: Function used to invoke executables on the system.
        For example 'mpirun', 'mpiexec', 'srun', 'ibrun'
    :type ntask_max: int
    :param ntask_max: limit the number of concurrent tasks in a given array job
    :type tasktime: float
    :param tasktime: maximum job time in minutes for each job spawned by
        the SeisFlows master job during a workflow. These include, e.g.,
        running the forward solver, adjoint solver, smoother, kernel combiner.
        All spawned tasks receive the same task time. Fractions of minutes
        acceptable.
    :type environs: str
    :param environs: Optional environment variables to be provided in the
        following format VAR1=var1,VAR2=var2... Will be set using
        os.environs

    Paths
    -----
    :type path_container: str
    :param path_container: path to the Docker Image that contains adjTomo
        software package
    ***
    """
    __doc__ = Workstation.__doc__ + __doc__

    def __init__(self, title=None, mpiexec="", ntask_max=None,
                 tasktime=1, environs="", singularity_exec="singularity",
                 path_container=None, **kwargs):
        """Instantiate the Cluster System class"""
        super().__init__(**kwargs)

        if title is None:
            self.title = os.path.basename(os.getcwd())
        else:
            self.title = title
        self.mpiexec = mpiexec
        self.ntask_max = ntask_max or nproc() - 1  # -1 because master job
        self.tasktime = tasktime
        self.environs = environs or ""
        self.singularity_exec = singularity_exec

        self.path["container"] = path_container
        self.path["run"] = os.path.join(self.path.scratch, "run")

    # def check(self):
    #     """
    #     Checks path and parameter validity
    #     """
    #     super().check()
    #
    #     assert(self.path.container and os.path.exists(self.path.container)), (
    #         f"`Singularity` System class requires a Docker Image specified in "
    #         f" path `path_container`"
    #     )
    #
    #     # # Check that 'singularity' exists on system
    #     try:
    #         subprocess.run(f"which {self.singularity_exec}", check=True,
    #                        shell=True, stdout=subprocess.DEVNULL)
    #     except subprocess.CalledProcessError:
    #         logger.critical(f"Singularity exectuable {self.singularity_exec} "
    #                         f"not found on sytem (using `which`), but is "
    #                         f"required to run a `Singularity` System. Please "
    #                         f"check your parameter: `singularity_exec`.")
    #         sys.exit(-1)

    def setup(self):
        """
        Copies 'submit' and 'run' .py scripts from the repository into the
        working directory so that the User can run these scripts directly.
        This is a manual step in order to allow Users to run with a container
        without using native environment commands (e.g., sbatch) from inside a
        container.
        """
        super().setup()

        # Location of scripts INSIDE the repository
        run_script = os.path.join(ROOT_DIR, "system", "runscripts", "run")
        cp(run_script, self.path.run)

    @property
    def run_call_header(self):
        """
        The run call defines the Singularity wrapper which executes run calls
        using the Docker image. It also binds the current working directory
        inside the container so that we can write back to the local filesystem.

        .. note::
            Generalized `cluster` returns empty string but child system
            classes will need to overwrite the submit call.

        :rtype: str
        :return: the system-dependent portion of a run call
        """
        return (f"{self.singularity_exec} exec -c "
                f"--bind $(pwd):/work {self.path.container} "
                f"bash -c 'cd /work;")

    def submit(self, workdir=None, parameter_file="parameters.yaml"):
        """
        Submits the main workflow job as a serial job submitted directly to
        the system that is running the master job

        :type workdir: str
        :param workdir: path to the current working directory
        :type parameter_file: str
        :param parameter_file: parameter file name used to instantiate the
            SeisFlows package
        """
        workflow = import_seisflows(workdir=workdir or self.path.workdir,
                                    parameter_file=parameter_file)
        workflow.check()
        workflow.setup()
        workflow.run()

    def run(self, funcs, single=False, **kwargs):
        """
        Runs tasks multiple times in parallel by submitting NTASK new jobs to
        system. The list of functions and its kwargs are saved as pickles files,
        and then re-loaded by each submitted process with specific environment
        variables. Each spawned process will run the list of functions.

        :type funcs: list of methods
        :param funcs: a list of functions that should be run in order. All
            kwargs passed to run() will be passed into the functions.
        :type single: bool
        :param single: run a single-process, non-parallel task, such as
            smoothing the gradient, which only needs to be run by once.
            This will change how the job array and the number of tasks is
            defined, such that the job is submitted as a single-core job to
            the system.
        :type run_call: str
        :param run_call: the call used to submit the run script. If None,
            attempts default run call which should be suited for the given
            system. Can be overwritten by child classes to involve other
            arguments
        """
        # Single tasks only need to be run one time, as `TASK_ID` === 0
        ntasks = {True: 1, False: self.ntask}[single]
        funcs_fid, kwargs_fid = pickle_function_list(functions=funcs,
                                                     path=self.path.scratch,
                                                     **kwargs)
        logger.info(f"running functions {[_.__name__ for _ in funcs]} on "
                    f"system {self.ntask} times")

        # Create the run call which will simply call an external Python script
        # e.g., run --funcs func.p --kwargs kwargs.p --environment ...
        run_call = " ".join([
            f"{self.run_call_header}",
            f"{self.path.run}",
            f"--funcs {funcs_fid}",
            f"--kwargs {kwargs_fid}",
            f"--environment {self.environs}"
            f"'"  # <- important, closes the bash command started in header
        ])
        logger.debug(run_call)

        # !!! This is SLURM dependent. Will need to change if other scheduler
        if single:
            modify_run_call_single_proc(run_call)

        # Write the run call inside a script to make it cleaner for User to call
        runscript = os.path.join(self.path.scratch, "run_task.sh")
        with open(runscript, "w") as f:
            f.write("#!/bin/bash -e\n")  # shebang
            f.write(run_call)

        input(msg.cli(f"sh {os.path.relpath(runscript)}",
                      header="task run call", border="=",
                      items=["> please execute the above command",
                             "> hit `enter` when all jobs are complete",
                             "> monitor jobs using `sacct` or `squeue`"])
              )
