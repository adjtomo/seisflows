#!/usr/bin/env python3
"""
Shadowfax is a University of Alaska Fairbanks (UAF) high performance computer,
operated by the Geophysical Institute's Research Computing Systems (RCS).
Shadowfax operates on a SLURM workload manager and therefore overloads the SLURM
System module.

.. notes::

    - Shadowfax has 96 cores split onto 2 sockets with 24 cores per socket
    - Also has 4x Nvidia A40 GPUs 
    - Running RHEL Fedora 8.7 (Oopta)
"""
import os
from datetime import timedelta

from seisflows.system.slurm import Slurm


class Shadowfax(Slurm):
    """
    System Shadowfax
    --------------
    University of Alaska Fairbanks HPC ShadowFax, SLURM based system

    Parameters
    ----------
    :type partition: str
    :param partition: Chinook has various partitions which each have their
        own number of cores per compute node. Available are: analysis, t1small,
        t2small, t1standard, t2standard, gpu
    :type gpu: int
    :param gpu: if set, sets the number of GPUs to use for simulation work and
        ignores the parameter `ntask`. If not set, defaults to CPU mode whichi
        requests `ntask` CPUs per simulation

    :type submit_to: str
    :param submit_to: (Optional) partition to submit the main/master job which 
        is a serial Python task that controls the workflow. 

    Paths
    -----

    ***
    """
    __doc__ = Slurm.__doc__ + __doc__


    def __init__(self, mpiexec="mpiexec", partition="defq", gpu=False, 
                 submit_to=None, **kwargs):
        """Chinook init"""
        super().__init__(**kwargs)

        self.mpiexec = mpiexec
        self.partition = partition
        self.submit_to = submit_to or self.partition
        self.gpu = gpu

        self._partitions = {"defq": 24,}

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
            f"{self.slurm_args or ''}",
            f"--job-name={self.title}",
            f"--output={self.path.output_log}",
            f"--error={self.path.output_log}",
            f"--ntasks=1",
            f"--partition={self.submit_to}",
            f"--time={self._walltime}"
        ])
        return _call

    def run_call(self, executable="", single=False, array=None, tasktime=None):
        """
        The run call defines the SBATCH call which is used to run tasks during
        an executing workflow. Like the submit call its arguments are dictated
        by the given system. Run calls are modified and called by the `run`
        function

        :type executable: str
        :param exectuable: the actual exectuable to run within the SBATCH 
            directive. Something like './script.py'
        :type array: str
        :param array: overwrite the `array` variable to run specific jobs. If
            not provided, then we will run jobs 0-{ntask}%{ntask_max}. Jobs 
            should be submitted in the format of a SLURM array string, 
            something like: 0,1,3,5 or 2-4,8-22
        :type single: bool
        :param single: flag to get a run call that is meant to be run on the
            mainsolver (ntask==1), or run for all jobs (ntask times). Examples
            of single process runs include smoothing, and kernel combination
        :rtype: str
        :return: the system-dependent portion of a run call
        """
        array = array or self.task_ids(single=single)  # get job array str
        if tasktime is None:
            # 0 is there just for initialization. `tasktime` will superceded
            tasktime = str(timedelta(minutes=0 or self.tasktime))
        else:
            tasktime = str(timedelta(minutes=tasktime))

        # Determine if this is a single-process or array job
        if single:
            ntasks = 1
            env = "SEISFLOWS_TASKID=0," 
        else:
            ntasks = self.nproc
            env = ""
       
        _call = " ".join([
             f"sbatch",
             f"{self.slurm_args or ''}",
             f"--job-name={self.title}",
             f"--ntasks={ntasks:d}",
             f"--partition={self.partition}",
             f"--time={tasktime}",
             f"--output={os.path.join(self.path.log_files, '%A_%a')}",
             f"--array={array}",
             f"--parsable",
             f"{executable}",  # <-- The actual script/program to run goes here
             f"--environment {env}{self.environs or ''}"
        ])
        # Overwrite if using GPU mode
        if self.gpu:
            _call.replace(f"--ntasks={ntasks:d}", f"--gres=gpu:{self.gpu}")

        return _call
