#!/usr/bin/env python3
"""
Frontera is one of the Texas Advanced Computing Center (TACC) HPCs.
https://frontera-portal.tacc.utexas.edu/

.. note:: Frontera Caveat 1  
    On TACC Systems you cannot submit 'sbatch' from compute nodes. Work around:
    SSHs from compute node to login node, activate conda environemtn, submit 
    sbatch script. This requires knowing the User name, conda environment name,
    and ensuring SSH keys are available. Thanks to Ian Wang for the suggestion
    to SSH around the problem.

    Essentially we are running, from the compute node:
    $ ssh user@hostname 'conda activate env; sbatch --arg=val run_function.sh'

.. note:: Frontera Caveat 2
    TACC does not allow the '--array' option, which SeisFlows uses to submit
    multiple jobs in a single SBATCH command. To work around this, the Frontera
    module submits jobs one by one.
"""
import os
import sys
import subprocess
from seisflows import ROOT_DIR, logger
from seisflows.tools import msg
from seisflows.tools.config import pickle_function_list
from seisflows.system.slurm import (Slurm, query_job_states, BAD_STATES, 
                                    check_job_status_list)


class Frontera(Slurm):
    """
    System Frontera
    --------------
    Texas Advanced Computing Center HPC Frontera, SLURM based system

    Parameters
    ----------
    :type user: str
    :param user: User's username on TACC systems. Can be determined by 'whoami'
        or will be gathered from the 'USER' environment variable. Used for
        internal ssh'ing from compute nodes to login nodes.
    :type conda_env: str
    :param conda_env: name of the Conda environment in which you are running 
        SeisFlows. Defaults to environment variable 'CONDA_DEFAULT_ENV'. Used
        to activate the conda environment AFTER ssh'ing from compute to login
        node, to ensure that the newly submitted job has access to the SeisFlows
        environment
    :type partition: str
    :param partition: Chinook has various partitions which each have their
        own number of cores per compute node. Available are: small, normal,
        large, development, flex
    :type allocation: str
    :param allocation: Name of allocation/project on the Frontera system.
        Required if you have more than one active allocation.

    Paths
    -----

    ***
    """
    def __init__(self, user=None, conda_env=None, partition="development", 
                 allocation=None, **kwargs):
        """Frontera init"""
        super().__init__(**kwargs)

        self.user = user or os.environ["USER"]  # alt. getpass.getuser()
        self.conda_env = conda_env or os.environ["CONDA_DEFAULT_ENV"]
        self.partition = partition
        self.allocation = allocation
        self.mpiexec = "ibrun"

        # See 'Frontera Caveat 1' note for why we need these calls
        self._ssh_call = f"ssh {self.user}@frontera.tacc.utexas.edu"
        self._conda_activate = f"conda activate {self.conda_env}"

        # Internally used check parameters. Because 'development' and 'large'
        # partitions do not allow >1 job per user, we cannot use them
        self._acceptable_partitions = ["small", "normal", "flex"]
        self._partitions = {"small": 28, "normal": 28, "large": 28, "flex": 28,
                            "development": 28}
        self._max_jobs = {"small": 20, "large": 1, "normal": 100, "flex": 15,
                          "development": 1}

        # Hard set `ntask_max` based on TACCS 'QOSMaxJobsPerUserLimit' 
        self.ntask_max = self._max_jobs[self.partition]

    def check(self):
        """
        Checks parameters and paths
        """
        super().check()

        assert(self.partition in self._acceptable_partitions), \
            f"Frontera `partition` must be in {self._acceptable_partitions}"

        assert(self._max_jobs[self.partition] > 1), (
            f"Frontera partition '{self.partition}' does not allow more than 1 "
            f"simultaneously running job, meaning SeisFlows will not work. "
            f"please choose a different partition"
            )

        if self.tasktime > 60 and self.partition == "flex":
            logger.warning("Frontera's 'Flex' partition may cancel jobs that "
                           "exceed 60 minute wall time. Consider choosing a "
                           "different partition if this may be a problem")

    @property
    def submit_call_header(self):
        """
        The submit call defines the SBATCH header which is used to submit a
        workflow task list to the system. It is usually dictated by the
        system's required parameters, such as account names and partitions.
        Submit calls are modified and called by the `submit` function.

        :rtype: str
        :return: the system-dependent portion of a submit call
        """
        _call = " ".join([
            f"sbatch",
            f"--job-name={self.title}",  # -J
            f"--partition={self.partition}",  # -p
            f"--output={self.path.output_log}",  # -o
            f"--error={self.path.output_log}",
            f"--nodes=1",  # -N
            f"--ntasks=1",  # -n
            f"--time={self._walltime}"  # -t
        ])
        if self.allocation is not None:
            _call = f"{_call} --allocation={self.allocation}"
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
            f"--job-name={self.title}",
            f"--partition={self.partition}",
            f"--output={os.path.join(self.path.log_files, '%A')}",
            f"--ntasks={self.nproc:d}",
            f"--nodes={self.nodes}",
            f"--time={self._tasktime}",
            f"--parsable"
        ])
        if self.allocation is not None:
            _call = f"{_call} --allocation={self.allocation}"
        return _call

    @staticmethod
    def _stdout_to_job_id(stdout):
        """                                                                      
        The stdout message after an SBATCH job is submitted. On Frontera, the
        standard message is preceded by a log message which looks like:

        ```
        -----------------------------------------------------------------
                   Welcome to the Frontera Supercomputer                 
        -----------------------------------------------------------------

        No reservation for this job
        --> Verifying valid submit host (login3)...OK
        --> Verifying valid jobname...OK
        --> Verifying valid ssh keys...OK
        --> Verifying access to desired queue (development)...OK
        --> Checking available allocation (EAR21042)...OK
        --> Verifying that quota for filesystem ... is at 3.87% allocated...OK
        --> Verifying that quota for filesystem ... is at  0.91% allocated...OK
        4738284

        ```
        :type stdout: str                                                        
        :param stdout: standard SBATCH response after submitting a job with the
            '--parsable' flag
        :rtype: str                                                              
        :return: a matching job ID. We convert str->int->str to ensure that      
            the job id is an integer value (which it must be)     
        """
        job_id = stdout.split("OK")[-1].strip().split(";")[0]
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
        Runs task multiple times in embarrassingly parallel fasion on Frontera.
        Executes the list of functions (`funcs`) NTASK times with each  
        task occupying NPROC cores.                                              
                                                                                 
        .. note::                                                                
            Completely overwrites the `Slurm.run()` command                    

         TODO 
            * can we ssh once or do we have to do it for each process?
            * the ssh command prints the ssh prompt to the log file, how do
              we supress that?
                                                                                 
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
                                                                                 
        # Run call is slightly different for Frontera, which needs ssh and 
        # cannot run array jobs
        job_ids = []
        for taskid in range(_ntask):
            run_call = " ".join([
                f"{self.run_call_header}",
                f"{os.path.join(ROOT_DIR, 'system', 'runscripts', 'run')}",
                f"--funcs {funcs_fid}",
                f"--kwargs {kwargs_fid}",
                f"--environment {self.environs or ''},SEISFLOWS_TASKID={taskid}"
            ])                    
            # Need to wrap run call in quotes so it can go through ssh
            run_call = f"{self._ssh_call} '{self._conda_activate}; {run_call}'"
            if taskid == 0:
                logger.debug(run_call)
                                                                                 
            stdout = subprocess.run(run_call, stdout=subprocess.PIPE,
                                    text=True, shell=True).stdout
            job_ids.append(self._stdout_to_job_id(stdout))

        # Monitor the job queue until all jobs have completed, or any one fails  
        try:                                                                     
            status = check_job_status_list(job_ids)
        except FileNotFoundError:                                                
            logger.critical(f"cannot access job information through 'sacct', "   
                            f"waited 50s with no return, please check job "      
                            f"scheduler and log messages")                       
            sys.exit(-1)                                                         
                                                                                 
        if status == -1:  # Failed job                                           
            logger.critical(                                                     
                msg.cli(f"Stopping workflow. Please check logs for details.",    
                        items=[f"TASKS:   {[_.__name__ for _ in funcs]}",        
                               f"SBATCH:  {run_call}"],                          
                        header="slurm run error", border="=")    
            )                                                                    
            sys.exit(-1)                                                         
        else:                                                                    
            logger.info(f"{self.ntask} tasks finished successfully") 


