#!/usr/bin/env python3
"""
Frontera is one of the Texas Advanced Computing Center (TACC) HPCs.
https://frontera-portal.tacc.utexas.edu/

.. note::
    One caveat of the TACC Systems is that you cannot submit 'sbatch' from 
    compute nodes, which is how SeisFlows operates. To work around this, 
    the run call SSHs from the compute node to the login node to submit the
    sbatch script. This requires knowing the User name, and that SSH keys
    are available. Thanks to Ian Wang for the suggestion.
"""
import os
import sys
from seisflows import logger
from seisflows.system.slurm import Slurm


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
    def __init__(self, user=None, partition="development", allocation=None, 
                 **kwargs):
        """Frontera init"""
        super().__init__(**kwargs)

        self.user = user or os.environ["USER"]  # alt. getpass.getuser()
        self.partition = partition
        self.allocation = allocation
        self.mpiexec = "ibrun"

        # See note in file docstring for why we need this SSH call
        self._ssh_call = f"ssh {self.user}@frontera.tacc.utexas.edu"
        self._partitions = {"small": 28, "normal": 28, "large": 28,
                            "development": 28, "flex": 28}


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
            f"{self._ssh_call}",
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
            f"{self._ssh_call}",
            f"sbatch",
            f"{self.slurm_args or ''}",
            f"--job-name={self.title}",
            f"--partition={self.partition}",
            f"--output={os.path.join(self.path.log_files, '%A_%a')}",
            f"--ntasks={self.nproc:d}",
            f"--nodes={self.nodes}",
            f"--array=0-{self.ntask - 1 % self.ntask_max}",
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
