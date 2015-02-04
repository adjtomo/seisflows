
Overview
========

SeisFlows is a Python adjoint tomography and full waveform inversion package designed to be flexible enough to use in a variety of contexts, including scientific research.

To provide this flexibility, users are offered choices in each of the following categories: workflow, system, solver, optimization, preprocessing, and postprocessing.  Within each category, different classes are interchangeable. If desired functionality is missing from the main package, users can customize default classes by overloading methods, or contribute their own classes.  This combination of extensibility, modular design, and object oriented programming allows multiple users to work productively within the same framework.

To illustrate how it works, consider an example from regional tomography.  A typical regional setup involves a 3D Cartesian solver run on a PBS cluster.  Under SeisFlows, if the study area expands, users can replace the 3D Cartesian solver with a 3D spherical solver.  If a faster SLURM cluster comes online, users can substitute the PBS system interface for a SLURM system interface. If a new type of data becomes available, users can modify the misfit function by overloading appropriate methods.  


Installation
============

To install Seisflows, first clone the repository::

    git clone github.com/PrincetonUniversity/seisflows


Then set environment variables. If using bash, add the following lines to ``.bash_profile`` or ``.bashrc``::

    export PATH=$PATH:/path/to/seisflows/scripts
    export PYTHONPATH=$PYTHONPATH:/path/to/seisflows
    export PYTHONUNBUFFERED="TRUE"


Software Prerequisites
----------------------

SeisFlows requires Python 2.7, NumPy >1.6, and SciPy >0.12. Forward modeling software is also a prerequisite; see :ref:`solver` for more information.


Hardware Prerequisites
----------------------

Access to a computer cluster is required for most applications.  Base classes are provided for several common cluster configurations, including PBS and SLURM.  Nonstandard configurations can often be accommodated through modifications to one of the base classes; see :ref:`system` for details.


.. _job_submission:

Job Submission
==============

Each job must be submitted from a `working directory`.  Within a working directory, users must supply two input files, ``paths.py`` and ``parameters.py``. Output files, by default, are written to the working directory, along with scratch files created by the solver and optimization routines. Different output and scratch directories can be specified by adding or modifying entries in ``paths.py``.

``parameters.py`` contains a list of parameter names and values. Prior to a job being submitted, parameters are checked so that errors can be detected without loss of queue time or wall time. Parameters are stored in a dictionary that is accessible from anywhere in the Python code. By convention, all parameter names must be upper case. Parameter values can be floats, integers, strings or any other Python data type. Parameters can be listed in any order.

``paths.py`` contains a list of path names and values. Prior to a job being submitted, paths are checked so that errors can be detected without loss of queue time or wall time. Paths are stored in a dictionary that is accessible from anywhere in the Python code. By convention, all names must be upper case, and all values must be absolute paths. Paths can be listed in any order.

Once a working directory and input files have been created, users can type ``sfrun`` from within the working directory to submit a job. If the ``serial`` system configuration is specified in ``parameters.py``, the job will begin executing immediately. If ``pbs`` or ``slurm`` configurations are specified, the job will run when resources become available. Once the job starts running, status information will be displayed either to the terminal or to the file ``output.log``.



.. _solver:

Solver Configuration
====================

SeisFlows includes Python interfaces for SPECFEM2D, SPECFEM3D, and SPECFEM3D_GLOBE.  While the Python interfaces are part of the SeisFlows package, the solver source code itself is located in other repositories and must be downloaded separately.  

After downloading the solver source code, users must configure and compile it, following the instructions in the solver user manual. Summarized briefly, the configuration and compilation procedure is as follows:

Prior to compilation, users need to run the ``configure`` script and prepare input files such as

- parameter file

- source file

- stations file.

The result of compilation is a set of binaries files, including

- mesher binary

- solver binary

- smoothing utility

- summing utility.


Note that if input files change, binary files may need to be recompiled.

After compilation, input files must be gathered together in one directory and binary files in another.  The absolute paths to input file and binary file directories must be given in ``paths.py`` as follows::

    SOLVER_INPUT = '/path/to/solver/input/files'
    SOLVER_BINARIES = '/path/to/solver/binary/files'


Solver Integration
------------------

Integration of the solver with the other workflow components can be challenging. Here we try to give an idea of the issues involved from both a developer and a user standpoint.

- Solver computations account for most of the cost of an inversion. As a result, the solver must be written in an efficient compiled language, and wrappers must be written to integrate the compiled code with other software components. 

- Because compilation of SPECFEM2D, SPECFEM3D, or SPECFEM3D_GLOBE depends on input files, there is currently no mechanism for automatically compiling binary files. Users must prepare their own SPECFEM input files and then follow the procedure from the SPECFEM documentation to compile binary files.

- As described :ref:`above <job_submission>`, SeisFlows uses its own unique input files to determine runtime settings.  Problems could arise if parameters from SeisFlows input files conflict with parameters from SPECFEM input files. Users must make sure that there are no conflicts between SeisFlows parameters and solver parameters.

- In the solver routines, it is natural to represent velocity models as dictionaries, with different keys corresponding to different material parameters.  In the optimization routines, it natural to represent velocity models as vectors. To convert back and forth between these two representations, a pair of utility functions--``split`` and ``merge``--are included in each in the SeisFlows package as part of the solver interfaces.


Writing Custom Solver Interfaces
--------------------------------

Besides SPECFEM2D, SPECFEM3D, and SPECFEM3D_GLOBE, SeisFlows can interface with other solvers capable of running forward and adjoint simulations. For information about writing custom solver interfaces, see :ref:`developer`.


.. _system:

System Configuration
====================

SeisFlows can run on SLURM, PBS TORQUE, and PBS Pro clusters.  For debugging, an option to run simulations in serial is also provided.  

While there are many similarities between job management systems, there are also many differences.  Our approach to such differences is to try to hide them behind a consistent Python interface.  For example, by creating a thin Python layer over system commands such as ``qsub`` on PBS or ``sbatch`` on SLURM, it is possible to abstract the machinery for submitting and managing jobs.

Besides different job submission and management systems, different filesystem configurations may exist as well.  Filesystem settings can be adjusted by modifying values in the ``PATH`` dictionary, which is populated from ``paths.py``.  Output files and temporary files, by default, are written to the working directory.  If a value for ``PATH.GLOBAL`` is supplied, temporary files are written there instead.  If each compute node has its own local filesystem and if a value for ``PATH.LOCAL`` is supplied, some temporary files will be written to ``PATH.LOCAL`` and others to ``PATH.GLOBAL``.

As the size of an inversion grows, scalability and fault tolerance become increasingly important.  If a single forward simulation spans more than one node, users must select ``pbs_lg_job`` or ``slurm_lg_job`` system configurations in ``parameters.py``.  If a forward simulation fits onto a single node, users must select ``pbs_lg_job`` or ``slurm_sm_job`` instead.


Heavyweight Solutions
---------------------

In writing system interfaces, the approach taken by SeisFlows developers has been to write lightweight Python wrappers on top of PBS and SLURM commands.  For some cases involving nonstandard cluster configurations or restrictive usage policies, heavyweight solutions may be required instead.  Users are referred to distributed computing projects such as SAGA or PATHOS for ideas.



.. _developer:

Developer Reference
===================

To allow classes to work with one another, each class must conform to an established interface.  In practice, this means each class must implement specified methods, listed below, with specified input and output.

``solver`` classes must implement

- check

- setup

- eval_func

- eval_grad

- forward

- adjoint

- load

- save

- split

- merge


``system`` classes must  implement

- check

- submit

- run


``preprocess`` classes must implement

- check

- setup

- prepare_eval_grad

- process_traces

- write_residuals


``postprocess`` classes must implement

- check

- setup

- process_kernels


``optimize`` classes must implement

- check

- setup

- compute_direction

- compute_step

- initialize_search

- finalize_search

- search_status


``workflow`` classes must implement

- check

- main


In the above list, ``setup`` methods are generic methods, called from the ``main`` workflow script and meant to provide users the flexibility to perform any required setup tasks. ``check`` methods are the default mechanism for parameter declaration and checking and are called just once, prior to a job being submitted through the scheduler.

Besides required methods, classes may include any number of private methods or utility functions.


System Interfaces
-----------------

A list of available system interface classes follows. By hiding environment details behind a python interface layer, these classes provide a consistent command set across different computing environments.

PBS_SM - For small inversions on PBS clusters. All resources are allocated at the beginning and all simulations are run at the same time, within a single job. Because of limitations of pbsdsh, individual wavefield simulations cannot span more than one core.

PBS_LG - For large inversions on PBS clusters. The work of the inversion is divided between multiple jobs, which are coordinated by a single long-running master job. Resources are allocated on a per simulation basis.

SLURM_SM - For small inversions on SLURM clusters. All resources are allocated at the beginning and all simulations are run at the same time, within a single job. Individual wavefield simulations can span more than one core, but span more than one node.

SLURM_LG - For large inversions on SLURM clusters. The work of the inversion is divided between multiple jobs, which are coordinated by a single long-running master job. Resources are allocated on a per simulation basis.

SLURM_XL - For very large inversions on SLURM clusters. In addition to the features of SLURM_LG, provides fault tolerence: Tasks that end in failure or timeout are automatically resumbitted. (Can be dangerous to use on code that is not well tested.)

SERIAL - Tasks that are normally run in parallel, such as multiple wavefield simulations, are carried out one at a time. Useful for debugging, among other things.

