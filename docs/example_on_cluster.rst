Running Examples on a Cluster
=============================

Here we detail how a User might transition from developing a 2D example
problem on their workstation, to performing large-scale inversion on a
cluster. In this notebook we show an example running on the New Zealand
eScience Infrastructure HPC, named Maui, but is meant to provide a
generalizable approach for running SeisFlows on clusters.

Example Setup
-------------

We first set up our working directory using the example setup shown in the `SPECFEM2D example page <specfem2d_example.html>`__. This ensures that we have our initial and final models, and a properly set parameter file that can be used for our inversion.

.. code:: ipython3

    # This is an empty working directory
    %cd /home/bchow/Work/scratch 


.. parsed-literal::

    /home/bchow/Work/scratch


.. code:: bash

    ! ln -s /home/bchow/REPOSITORIES/specfem2d .  # place SPECFEM2D repository in the working directory
    ! seisflows examples setup 2  # run example setup but do not `submit` workflow

.. code:: ipython3

    ! ls 


.. parsed-literal::

    parameters.yaml  specfem2d  specfem2d_workdir


Module ‘swap’
-------------

As we saw in the example tutorial, the ``System`` module for this
example problem is set as ‘Workstation’, which is meant to run the
workflow in serial directly on the system that submits it. For clusters
this means we would run our entire inversion on the login node.

.. code:: ipython3

    ! seisflows par system


.. parsed-literal::

    system: workstation


To ‘swap’ out the ``System`` module for a cluster-specific class, we can
use the ``seisflows swap`` command, which replaces one module for
another without affecting the other modules. This is very helpful if you
have a completed parameter file and do not want to copy-paste all the
edited parameter just to change out a module. The rubric for running
``seisflows swap`` can be found in the help message:

.. code:: ipython3

    ! seisflows swap -h


.. parsed-literal::

    usage: seisflows swap [-h] [module] [classname]
    
    During workflow development, it may be necessary to swap between different
    sub-modules (e.g., system.workstation -> system.cluster). However this would
    typically involving re-generating and re-filling a parameter file. The 'swap'
    function makes it easier to swap parameters between modules.
    
    positional arguments:
      module      Module name to swap
      classname   Classname to swap to
    
    optional arguments:
      -h, --help  show this help message and exit


You can check available names by running ``seisflows print modules``.
Here we want to swap out our ``System`` module from ‘Workstation’ to
‘Maui’, which defines how SeisFlows interacts with the SLURM-based
system, Maui.

.. code:: ipython3

    ! seisflows print modules


.. parsed-literal::

                                   SEISFLOWS MODULES                                
                                   /////////////////                                
    '-': module, '*': class
    
    - workflow
        * forward
        * inversion
        * migration
        * test_flow
    - system
        * chinook
        * cluster
        * frontera
        * lsf
        * maui
        * slurm
        * workstation
    - solver
        * specfem
        * specfem2d
        * specfem3d
        * specfem3d_globe
    - preprocess
        * default
        * pyaflowa
    - optimize
        * LBFGS
        * NLCG
        * gradient


.. code:: ipython3

    ! seisflows swap system maui


.. parsed-literal::

    L-BFGS optimization requires 'backtrack'ing line search. Overwriting 'bracket'


We can see now that the parameter file has swapped out the ‘Workstation’
System module for the ‘Maui’ System module, which contains its own set
of parameters that must be filled out by the User.

.. code:: ipython3

    ! head -235 parameters.yaml | tail -n 110 


.. parsed-literal::

    # =============================================================================
    #
    #    Workstation System
    #    ------------------
    #    Defines foundational structure for System module. When used standalone, 
    #    runs tasks in serial on a local machine.
    #
    #    Parameters
    #    ----------
    #    :type ntask: int
    #    :param ntask: number of individual tasks/events to run during workflow.
    #        Must be <= the number of source files in `path_specfem_data`
    #    :type nproc: int
    #    :param nproc: number of processors to use for each simulation
    #    :type log_level: str
    #    :param log_level: logger level to pass to logging module.
    #        Available: 'debug', 'info', 'warning', 'critical'
    #    :type verbose: bool
    #    :param verbose: if True, formats the log messages to include the file
    #        name, line number and message type. Useful for debugging but
    #        also very verbose.
    #
    #        
    #    Cluster System
    #    ------------------
    #    Generic or common HPC/cluster interfacing commands
    #
    #    Parameters
    #    ----------
    #    :type title: str
    #    :param title: The name used to submit jobs to the system, defaults
    #        to the name of the current working directory
    #    :type mpiexec: str
    #    :param mpiexec: Function used to invoke executables on the system.
    #        For example 'mpirun', 'mpiexec', 'srun', 'ibrun'
    #    :type ntask_max: int
    #    :param ntask_max: limit the number of concurrent tasks in a given array job
    #    :type walltime: float
    #    :param walltime: maximum job time in minutes for the master SeisFlows
    #        job submitted to cluster. Fractions of minutes acceptable.
    #    :type tasktime: float
    #    :param tasktime: maximum job time in minutes for each job spawned by
    #        the SeisFlows master job during a workflow. These include, e.g.,
    #        running the forward solver, adjoint solver, smoother, kernel combiner.
    #        All spawned tasks receive the same task time. Fractions of minutes
    #        acceptable.
    #    :type environs: str
    #    :param environs: Optional environment variables to be provided in the
    #        following format VAR1=var1,VAR2=var2... Will be set using
    #        os.environs
    #
    #        
    #    System Slurm
    #    ------------------
    #    Interface for submitting and monitoring jobs on HPC systems running the 
    #    Simple Linux Utility for Resource Management (SLURM) workload manager.
    #
    #    Parameters
    #    ----------
    #    :type slurm_args: str
    #    :param slurm_args: Any (optional) additional SLURM arguments that will
    #        be passed to the SBATCH scripts. Should be in the form:
    #        '--key1=value1 --key2=value2"
    #
    #        
    #    System Maui
    #    -----------
    #    New Zealand Maui-specfic modifications to base SLURM system
    #
    #    Parameters
    #    ----------
    #    :type account: str
    #    :param account: Maui account to submit jobs under, will be used for the
    #        '--account' sbatch argument
    #    :type cpus_per_task: int
    #    :param cpus_per_task: allow for multiple cpus per task, i.e,.
    #        multithreaded jobs
    #    :type cluster: str
    #    :param cluster: cluster to submit jobs to. Available are Maui and
    #        Mahuika
    #    :type partition: str
    #    :param partition: partition of the cluster to submit jobs to.
    #    :type ancil_cluster: str
    #    :param ancil_cluster: name of the ancilary cluster used for pre-
    #        post-processing tasks.
    #    :type ancil_partition: name of the partition of the ancilary cluster
    #    :type ancil_tasktime: int
    #    :param ancil_tasktime: Tasktime in minutes for pre and post-processing
    #        jobs submitted to Maui ancil.
    #
    #        
    # =============================================================================
    ntask: 1
    nproc: 1
    log_level: DEBUG
    verbose: False
    title: scratch
    mpiexec:  None
    ntask_max: 100
    walltime: 10
    tasktime: 1
    environs: SLURM_MEM_PER_CPU
    slurm_args:  None
    partition: nesi_research
    account: None
    cluster: maui
    cpus_per_task: 1
    ancil_cluster: maui_ancil
    ancil_partition: nesi_prepost
    ancil_tasktime: 1


’Check’ing parameter validity
-----------------------------

Most of the default values should be okay for our purposes, but it’s up
the User to read the docstrings and determine if any of the default
values should be changed. If we run ``seisflows check`` we can check if
any of our parameters are incorrectly set.

.. code:: ipython3

    ! seisflows check


.. parsed-literal::

    
    ================================================================================
                                    PARAMETER ERRROR                                
                                    ////////////////                                
    System 'Maui' requires parameter 'account'
    ================================================================================


The ``Maui`` System check function has told us that it requires that the
parameter ``account`` be set. Note that these requirements will change
between different clusters, which dictate different SLURM parameters
when submitting jobs. We can specify the account parameter using the
``seisflows par`` command.

.. code:: ipython3

    ! seisflows par account gns03247


.. parsed-literal::

    account: null -> gns03247


.. code:: ipython3

    ! seisflows check

The ``seisflows check`` function has passed and we have succesfully
swapped out our System module with the ``Maui`` child class. Under the
hood, this class should take care of all the required interactions
between SeisFlows and the compute node. Now all that is left to do is to
run ``seisflows submit``, which should submit the master job to the
system and run our inversion on compute nodes.

TestFlow: Live testing SeisFlows on System
------------------------------------------

While developing, debugging or testing SeisFlows on System, it is not
ideal to submit simulation-based workflows, as these eat large amounts
of computational resources and may introduce problems of there own.

Here we introduce ‘TestFlow’, a SeisFlows workflow that runs simple test
functions on a cluster. This allows Users to check if SeisFlows can
appropriately interact with the HPC system with tasks like submitting
jobs, monitoring the job queue and catching failing jobs.

Below we show how to set up TestFlow for our test bed HPC, Maui. First
we generate a template parameter file and set the modules appropriately.

.. code:: ipython3

    # This is an empty working directory
    %rm -r /home/bchow/Work/scratch 
    %mkdir /home/bchow/Work/scratch 
    %cd /home/bchow/Work/scratch 


.. parsed-literal::

    shell-init: error retrieving current directory: getcwd: cannot access parent directories: No such file or directory
    /home/bchow/Work/scratch


.. code:: ipython3

    # Generate a template parameter file
    ! seisflows setup -f


.. parsed-literal::

    creating parameter file: parameters.yaml


.. code:: ipython3

    # Set the modules appropriately
    ! seisflows par workflow test_flow
    ! seisflows par system maui  # we want to test SeisFlows on Maui
    ! seisflows par solver null  # currently test_flow does not test solver
    ! seisflows par preprocess null  # currently test_flow does not test preprocess
    ! seisflows par optimize null  # currently test_flow does not test optimize


.. parsed-literal::

    workflow: forward -> test_flow
    system: workstation -> maui
    solver: specfem2d -> null
    preprocess: default -> null
    optimize: gradient -> null


.. code:: ipython3

    # Dynamically fill out the parameter file
    ! seisflows configure

.. code:: ipython3

    ! head -48 parameters.yaml | tail -n 16


.. parsed-literal::

    # =============================================================================
    #
    #    TestFlow Workflow
    #    -------------
    #    Test individual sub-modules in a 'live' testing environment in order to
    #    ensure SeisFlows works appropriately given an established system and solver.
    #
    #    .. note::
    #        You do not need to set System parameters `ntask`, `nproc`, `tasktime`,
    #        `walltime`. These will be overwritten by the setup task.
    #
    #    Parameters
    #    ----------
    #
    #        
    # =============================================================================


As we can see above, the ``TestFlow`` workflow does not require any
input parameters, and will additionally automatically set some key
``System`` parameters to ensure that these tests are lightweight to
avoid long queue times. Under the hood the ``TestFlow`` workflow will:

1) Submit an array job to the system to test job submission capabilities
2) Submit a single job to the system which is intended to fail, this
   tests job queue monitoring as well as failed job catching.

Developers who are implementing new ``System`` classes (e.g., for new
clusters), can use TestFlow as foundation for their development and
debugging sessions. To run the ``TestFlow`` workflow you just need to
run ``seisflows submit``

.. code:: ipython3

    ! seisflows submit
