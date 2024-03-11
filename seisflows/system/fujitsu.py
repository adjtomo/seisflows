#!/usr/bin/env python3
"""
Fujitsu brand clusters run their own Job Management System (Fujitsu Technical 
Computing Suite). 

.. note::

    The nickname `PJM`, based on the batch job script directives, may be used 
    as a shorthand to refer to the Fujitsu job management system.

PJM is similar to PBS/TORQUE but has its own unique peculiarities. Systems like 
UToyko's Wisteria, and Fugaku (one of the top supercomputers in the world) run 
using the Fujitsu job scheduler.
"""
import os
import sys
import numpy as np
import time
import subprocess

from datetime import timedelta
from seisflows import logger
from seisflows.system.cluster import Cluster
from seisflows.tools import msg
from seisflows.tools.config import pickle_function_list, import_seisflows   


class Fujitsu(Cluster):
    """
    System Fujitsu
    --------------
    Interface for submitting and monitoring jobs on HPC systems running the 
    Fujitsu job management system, a.k.a PJM

    Parameters
    ----------
    :type pjm_args: str                                                        
    :param pjm_args: Any (optional) additional PJM arguments that will       
        be passed to the pjsub scripts. Should be in the form:                  
        '--key1=value1 --key2=value2"    

    Paths
    -----
    ***
    """
    __doc__ = Cluster.__doc__ + __doc__

    def __init__(self, ntask_max=100, pjm_args="", **kwargs):
        """
        Fujitsu-specific setup parameters

        :type ntask_max: int
        :param ntask_max: set the maximum number of simultaneously running 
            job processes that are submitted to a cluster at one time.
        """
        super().__init__(**kwargs)

        # Overwrite the existing 'mpiexec'
        if self.mpiexec is None:
            self.mpiexec = "mpiexec"
        self.ntask_max = ntask_max
        self.pjm_args = pjm_args

        # Must be filled in by child class
        self.group = None
        self.rscgrp = None
        self._rscgrps = {}

        # Convert walltime and tasktime to datetime str 'H:MM:SS'
        self._tasktime = str(timedelta(minutes=self.tasktime))
        self._walltime = str(timedelta(minutes=self.walltime))
    
        # Define PJM-dependent job states used for monitoring queue which
        # can be found in the User manual
        self._completed_states = ["END"]
        self._failed_states = ["CANCEL", "HOLD", "ERROR"]
        self._pending_states = ["QUEUED", "RUNNING"]

    def check(self):
        """
        Checks parameters and paths
        """
        super().check()

        assert(self.node_size is not None), (                                    
            f"Fujitsu system child classes require defining the node_size or "     
            f"the number of cores per node inherent to the compute system.")     
                                                                                 
        assert(self.rscgrp in self._rscgrps), \
            f"Cluster resource group must match {self._rscgrps}"  

    @property
    def nodes(self):
        """Defines the number of nodes which is derived from system node size"""
        _nodes = np.ceil(self.nproc / float(self.node_size))
        _nodes = _nodes.astype(int)
        return _nodes

    @property
    def node_size(self):
        """Defines the node size of a given cluster partition. This is a hard
        set number defined by the system architecture"""
        return self._rscgrps[self.rscgrp]

    @property
    def submit_call_header(self):
        """
        The submit call defines the PJM header which is used to submit a
        workflow task list to the system. It is usually dictated by the
        system's required parameters, such as account names and partitions.
        Submit calls are modified and called by the `submit` function.

        :rtype: str
        :return: the system-dependent portion of a submit call
        """
        _call = " ".join([
            f"pjsub",
            f"{self.pjm_args or ''}",
            f"-L rscgrp={self.rscgrp}",  # resource group
            f"-g {self.group}",  # project code
            f"-N {self.title}",  # job name
            f"-o {self.path.output_log}",  # write stdout to file
            f"-e {self.path.output_log}",  # write stderr to file
            f"-L elapse={self._walltime}",  # [[hour:]minute:]second
            f"-L node=1",
            f"--mpi proc=1",
        ])
        return _call

    @property
    def run_call_header(self):
        """
        The run call defines the PJM header which is used to run tasks during
        an executing workflow. Like the submit call its arguments are dictated
        by the given system. Run calls are modified and called by the `run`
        function

        :rtype: str
        :return: the system-dependent portion of a run call
        """
        _call = " ".join([
             f"pjsub",
             f"{self.pjm_args or ''}",
             f"-L rscgrp={self.rscgrp}",  # resource group
             f"-g {self.group}",  # project code
             f"-N {self.title}",  # job name
             f"-o {os.path.join(self.path.log_files, '%j')}", 
             f"-j",  # merge stderr with stdout
             f"-L elapse={self._tasktime}",  # [[hour:]minute:]second
             f"-L node={self.nodes}",
             f"--mpi proc={self.nproc}",
        ])
        return _call
    
    def submit(self, workdir=None, parameter_file="parameters.yaml", 
               direct=True):
        """
        Submit main workflow to the System. Two options are available,
        submitting a Python job directly to the system, or submitting a 
        subprocess.

        :type workdir: str                                                       
        :param workdir: path to the current working directory                    
        :type parameter_file: str                                                
        :param parameter_file: paramter file file name used to instantiate       
            the SeisFlows package    
        :type direct: bool
        :param direct: (used for overriding system modules) submits the main 
            workflow job directly to the login node as a Python process 
            (default). If False, submits the main job as a separate subprocess.
            Note that this is Fujitsu specific and main jobs should be run from
            interactive jobs run on compute nodes to avoid running jobs on
            shared login resources
        """
        if direct:
            workflow = import_seisflows(workdir=workdir or self.path.workdir,        
                                        parameter_file=parameter_file)               
            workflow.check()                                                         
            workflow.setup()                                                         
            workflow.run()    
        else:
            # e.g., submit -w ./ -p parameters.yaml                                  
            submit_call = " ".join([                                                 
                f"{self.submit_call_header}",
                f"{self.submit_workflow}",
            ])
                                                                                    
            logger.debug(submit_call)                                                
            try:                                                                     
                subprocess.run(submit_call, shell=True)                              
            except subprocess.CalledProcessError as e:                               
                logger.critical(f"SeisFlows master job has failed with: {e}")        
                sys.exit(-1)                                        

    @staticmethod                                                                
    def _stdout_to_job_id(stdout):                                               
        """                                                                      
        The stdout message after a PJSUB job is submitted, from which we get   
        the job number, example std output after submission:

        [INFO] PJM 0000 pjsub Job 1334958 submitted.

        :type stdout: str                                                        
        :param stdout: standard PJM output that is gathered from subprocess run
        :rtype: str                                                              
        :return: a matching job ID. We convert str->int->str to ensure that      
            the job id is an integer value (which it must be)                    
        :raises SystemExit: if the job id does not evaluate as an integer        
        """                                                                      
        job_id = str(stdout).split()[5].strip()                               
        try:                                                                     
            int(job_id)                                                          
        except ValueError:                                                       
            logger.critical(f"parsed job id '{job_id}' does not evaluate as an " 
                            f"integer, please check that function "              
                            f"`system._stdout_to_job_id()` is set correctly")    
            sys.exit(-1)                                                         
                                                                                 
        return job_id   

    def run(self, funcs, single=False, **kwargs):
        """
        Runs task multiple times in embarrassingly parallel fasion on a PJM
        cluster. Executes the list of functions (`funcs`) NTASK times with each
        task occupying NPROC cores.

        .. note::
        
            Completely overwrites the `Cluster.run()` command

        :type funcs: list of methods
        :param funcs: a list of functions that should be run in order. All
            kwargs passed to run() will be passed into the functions.
        :type single: bool
        :param single: run a single-process, non-parallel task, such as
            smoothing the gradient, which only needs to be run by once.
            This will change how the job array and the number of tasks is
            defined, such that the job is submitted as a single-core job to
            the system.
        """
        funcs_fid, kwargs_fid = pickle_function_list(funcs,
                                                     path=self.path.scratch,
                                                     **kwargs)
        if single:
            logger.info(f"running functions {[_.__name__ for _ in funcs]} on "
                        f"system 1 time")
            _ntask = 1
        else:
            logger.info(f"running functions {[_.__name__ for _ in funcs]} on "
                        f"system {self.ntask} times")
            _ntask = self.ntask

        # If no environs, ensure there is not trailing comma
        if self.environs:
            self.environs = f",{self.environs}"

        # Default Fujitsu command line input, can be overloaded by subclasses
        # Copy-paste this default run_call and adjust accordingly for subclass
        job_ids = []
        for taskid in range(_ntask):
            run_call = " ".join([
                f"{self.run_call_header}",
                # -x in 'pjsub' sets environment variables which are distributed
                # in the run script, see custom run scripts for example how
                # Ensure that these are comma-separated, not space-separated
                f"-x SEISFLOWS_FUNCS={funcs_fid},"  
                f"SEISFLOWS_KWARGS={kwargs_fid},"
                f"SEISFLOWS_TASKID={taskid}{self.environs},"
                f"GPU_MODE={int(bool(self.gpu))}",  # 0 if False, 1 if True
                f"{self.run_functions}",
            ])

            if taskid == 0:
                logger.debug(run_call)

            # Grab the job ids from each stdout message
            stdout = subprocess.run(run_call, stdout=subprocess.PIPE,
                                    text=True, shell=True).stdout
            job_ids.append(self._stdout_to_job_id(stdout))

        # Monitor the job queue until all jobs have completed, or any one fails
        try:
            status = self.monitor_job_status(job_ids)
        except FileNotFoundError:
            logger.critical(f"cannot access job information through 'pjstat', "
                            f"waited 50s with no return, please check job "
                            f"scheduler and log messages")
            sys.exit(-1)

        if status == -1:  # Failed job
            logger.critical(
                msg.cli(f"Stopping workflow. Please check logs for details.",
                        items=[f"TASKS:  {[_.__name__ for _ in funcs]}",
                               f"PJSUB:  {run_call}"],
                        header="PJM run error", border="=")
            )
            sys.exit(-1)
        else:
            logger.info(f"{self.ntask} tasks finished successfully")
            # Wait for all processes to finish and write to disk (if they do)
            # Moving on too quickly may result in required files not being avail
            time.sleep(5)

    def query_job_states(self, job_id, timeout_s=300, wait_time_s=30, 
                         _recheck=0):
        """
        Overwrites `system.cluster.Cluster.query_job_states`

        Queries completion status of array job by running the PJM cmd `pjstat`
        Because `pjstat` treats running and finished jobs separately, we need to
        query both the current running processes, and the history

        .. note::
            The actual command line call wil look something like this

            $ pjstat 
            Wisteria/BDEC-01 scheduled stop time: 2023/03/31(Fri) 09:00:00 (Remain: 28days 21:34:40)

            JOB_ID JOB_NAME STATUS PROJECT RSCGROUP START_DATE ELAPSE TOKEN NODEGPU
            1334861 STDIN END gr58interactive-o_n1 03/02 11:16:33  00:02:14 - 1 -

        .. note::
            `pjstat` flag options are described as follows:
            -H: get the history (non-running jobs)

        :type job_id: str
        :param job_id: main job id to query, returned from the subprocess.run 
            that ran the jobs
        :type timeout_s: float
        :param timeout_s: length of time [s] to wait for system to respond 
            nominally before rasing a TimeoutError. Sometimes it takes a while
            for job monitoring to start working.
        :type wait_time_s: float
        :param wait_time_s: initial wait time to allow System to initialize 
            jobs in the queue. I.e., how long do we  expect it to take before 
            the job we submitted shows up in the queue.
        :type _recheck: int
        :param _recheck: Used for recursive calling of the function. It can take 
            time for jobs to be initiated on a system, which may result in the 
            stdout of the 'sacct' command to be empty. In this case we wait and 
            call the function again. Rechecks are used to prevent endless loops 
            by putting a stop criteria
        :raises TimeoutError: if 'pjstat' does not return any output for ~1 min.
        """
        job_ids, job_states = [], []

        # Look at currently running jobs as well as completed jobs
        cmd = f"pjstat {job_id}"
        stdout = subprocess.run(cmd, stdout=subprocess.PIPE,
                                text=True, shell=True).stdout
        cmd = f"pjstat -H {job_id}"
        stdout += subprocess.run(cmd, stdout=subprocess.PIPE,
                                text=True, shell=True).stdout

        # Parse through stdout to get the job ID status
        for line in stdout.split("\n"):
            if line.startswith(str(job_id)):
                job_ids.append(line.split()[0])  # JOB_ID
                job_states.append(line.split()[2])  # STATUS 
        
        # Recursively re-check job state incase the job has not been 
        # instantiated in which cause 'stdout' is an empty string
        if not job_ids:
            _recheck += 1
            if _recheck > (timeout_s // wait_time_s):
                raise TimeoutError(f"cannot access job ID {job_id}")
            time.sleep(wait_time_s)
            self.query_job_states(job_id, _recheck=_recheck)

        return job_ids, job_states

