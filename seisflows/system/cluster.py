#!/usr/bin/env python3
"""
The Cluster class provides the core utilities interaction with HPC systems
which must be overloaded by subclasses for specific workload managers, or
specific clusters.
"""
import os
import dill
import subprocess
from seisflows import logger
from seisflows.config import ROOT_DIR
from seisflows.system.workstation import Workstation


class Cluster(Workstation):
    """
    [system.cluster] generic or common HPC/cluster interfacing commands

    :type title: str
    :param title: The name used to submit jobs to the system, defaults
        to the name of the current working directory
    :type mpiexec: str
    :param mpiexec: Function used to invoke executables on the system.
        For example 'mpirun', 'mpiexec', 'srun', 'ibrun'
    :type walltime: int
    :param walltime: maximum job time in minutes for the master SeisFlows
        job submitted to cluster
    :type tasktime: int
    :param tasktime: maximum job time in minutes for each job spawned by
        the SeisFlows master job during a workflow. These include, e.g.,
        running the forward solver, adjoint solver, smoother, kernel combiner.
        All spawned tasks receive the same task time.
    :type environs: str
    :param environs: Optional environment variables to be provided in the
        following format VAR1=var1,VAR2=var2... Will be set using
        os.environs
    """
    __doc__ = Workstation.__doc__ + __doc__

    def __init__(self, title=None, mpiexec="", walltime=10, tasktime=1,
                 environs="", **kwargs):
        """Instantiate the Cluster System class"""
        super().__init__(**kwargs)

        if title is None:
            self.title = os.path.basename(os.getcwd())
        else:
            self.title = title
        self.mpiexec = mpiexec
        self.walltime = walltime
        self.tasktime = tasktime
        self.environs = environs

    def _pickle_func_list(self, funcs, **kwargs):
        """
        Save a list of functions and their keyword arguments as pickle files.
        Return the names of the files for the run() function.

        .. note::
            The idea here is that we need this list of functions to be
            discoverable by a system separate to the one that defined them. To
            do this we can pickle Python objects on disk, and have the new
            system read in the pickle files and evaluate the objects. We use
            'dill' because Pickle can't serialize methods/functions

        :type funcs: list of methods
        :param funcs: a list of functions that should be run in order. All
            kwargs passed to run() will be passed into the functions.
        """
        # Save the instances that define the functions as a pickle object
        func_names = "_".join([_.__name__ for _ in funcs])  # unique identifier
        fid_funcs_pickle = os.path.join(self.path.scratch, f"{func_names}.p")

        with open(fid_funcs_pickle, "wb") as f:
            dill.dump(obj=funcs, file=f)

        # Save the kwargs as a separate pickle object
        fid_kwargs_pickle = os.path.join(self.path.scratch,
                                         f"{func_names}_kwargs.p")
        with open(fid_kwargs_pickle, "wb") as f:
            dill.dump(obj=kwargs, file=f)

        return fid_funcs_pickle, fid_kwargs_pickle

    def run(self, funcs, single=False, run_call=None, **kwargs):
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
            system
        """
        funcs_fid, kwargs_fid = self._pickle_func_list(funcs, **kwargs)
        logger.info(f"running functions {[_.__name__ for _ in funcs]} on "
                    f"system {self.ntask} times")

        if run_call is None:
            run_call = " ".join([
                f"{os.path.join(ROOT_DIR, 'system', 'runscripts', 'run')}",
                f"--funcs {funcs_fid}",
                f"--kwargs {kwargs_fid}",
                f"--environment SEISFLOWS_TASKID={{task_id}},{self.environs}"
            ])
        logger.debug(run_call)

        for task_id in range(self.ntask):
            logger.debug(f"running task id {task_id} "
                         f"(job {task_id + 1}/{self.ntask})")
            # Subprocess waits for the process to end before running the next
            subprocess.run(run_call.format(task_id=task_id), shell=True)
