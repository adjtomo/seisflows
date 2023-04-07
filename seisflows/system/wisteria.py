#!/usr/bin/env python3
"""
Wisteria is the University of Tokyo Fujitsu brand high performance computer.
Wisteria runs on the Fujitsu/PJM job scheduler.

.. note::

    - Wisteria has two node groups, Odyssey (compute nodes) and Aquarius 
      (data/learning nodes w/ GPU)
    - Odyssey has 7680 nodes with 48 cores/node
    - Aquarius has 45 nodes with 36 cores/node

.. note:: Wisteria Caveat 1     
                                                 
    On Wisteria you cannot submit batch jobs from compute nodes and you cannot
    SSH from compute nodes (Manual 5.13), so the master job must be 
    run from the login node or the pre-post node (Manual 5.2.3)

.. note:: Wisteria Caveat 2    
                                                  
    On Wisteria, the login node Conda environment is not inherited by compute
    nodes, so it requires custom `submit` and `run` script which first load
    the correct modules, and then run the corresponding script

.. note:: Wisteria Caveat 3

    On Wisteria, command line arguments for the `submit` and `run` script, 
    normally input like '--key value' interfere with the batch submission cmd
    `pjsub`. So instead we use the `pjsub` '-x' flag which allows us to set 
    environment variables. We use these in place of command line arguments
"""
import os
import subprocess
import sys
import time
from seisflows import ROOT_DIR, logger
from seisflows.tools.config import import_seisflows, pickle_function_list
from seisflows.system.fujitsu import Fujitsu, check_job_status_list


class Wisteria(Fujitsu):
    """
    System Wisteria
    ---------------
    University of Tokyo HPC Wisteria, running Fujitsu job scheduler

    Parameters
    ----------
    :type group: str
    :param group: User's group for allocating and charging resources. In the 
        pjsub script this is the '-g' option.
    :type rscgrp: str
    :param rscgrp: the resource group (i.e., partition) to submit jobs to. In 
        the pjsub script this is the '-L rscgrp' option.
        Available `rscgrp`s for Wisteria are:

        - debug-o: Odyssey debug, 30 min max, [1, 144] nodes available
        - short-o: Odyssey short, 8 hr. max, [1, 72] nodes available
        - regular-o: Odyssey regular, 24-48 hr. max, [1, 2304] nodes available
        - priority-o: Odyssey priority, 48 hr. max, [1, 288] nodes available

        - debug-a: Aquarius debug, 30 min max, [1, 1] nodes available
        - short-a: Aquarius short, 2 hr. max, [1, 2] nodes available
        - regular-a: Aquarius regular, 24-48 hr. max, [1, 8] nodes available

    Paths
    -----

    ***
    """
    __doc__ = Fujitsu.__doc__ + __doc__

    # Overwrites submit and run call locations to provide custom scripts that
    # are used to deal with the issue that the Conda environment is not 
    # automatically inherited on a compute node
    submit_workflow = os.path.join(ROOT_DIR, "system", "runscripts",
                                   "custom_submit-wisteria")   
    run_functions = os.path.join(ROOT_DIR, "system", "runscripts", 
                                 "custom_run-wisteria")   

    def __init__(self, user=None, group=None, rscgrp=None, **kwargs):
        """Wisteria init"""
        super().__init__(**kwargs)

        self.group = group
        self.rscgrp = rscgrp

        # Wisteria resource groups and their cores per node
        self._rscgrps = {
                "debug-o": 48, "short-o": 48, "regular-o": 48, "priority-o": 48,
                "debug-a": 48, "short-a": 48, "regular-a": 48
                }

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

        # Default Fujitsu command line input, can be overloaded by subclasses
        # Copy-paste this default run_call and adjust accordingly for subclass
        job_ids = []
        for taskid in range(_ntask):
            run_call = " ".join([
                f"{self.run_call_header}",
                # -x in 'pjsub' sets environment variables which are distributed
                # in the run script, see custom run scripts for example how
                f"-x SEISFLOWS_FUNCS={funcs_fid},SEISFLOWS_KWARGS={kwargs_fid},"
                f"SEISFLOWS_TASKID={taskid}",
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
            status = check_job_status_list(job_ids)
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

