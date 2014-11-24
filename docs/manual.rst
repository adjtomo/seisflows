
OVERVIEW
========

SeisFlows is an open source adjoint tomography and full waveform inversion package flexible enough to be used for a variety of purposes, including scientific research.

To provide this flexibility, users are offered choices in each of the following categories: workflow, system, solver, optimization, preprocessing, and postprocessing.  Within each category, different classes are interchangeable. If desired functionality is missing from the main package, users can customize default classes by overloading methods, or contribute their own classes.  This combination of modular design and object oriented programming style allows multiple users to work productively within the same framework.

To illustrate the how it all works, consider an example from regional tomography.  A typical set up might involve a 3D Cartesian solver and PBS cluster.  Under SeisFlows, if the study area expands, users can replace the 3D Cartesian solver with a 3D spherical solver.  If a faster computer cluster becomes available, users can replace the old system interface with a new one, say, SLURM.  All this is possible because of a system of simple, well defined interfaces between classes.



INSTALLATION
============

To install Seisflows, first clone the repository::

    git clone github.com/PrincetonUniversity/seisflows


Then set environment variables. If using bash, add the following lines to ``.bash_profile`` or ``.bashrc``::

    export PATH=$PATH:/path/to/seisflows/scripts
    export PYTHONPATH=$PYTHONPATH:/path/to/seisflows


Software Prerequisites
----------------------

SeisFlows requires Python 2.7, numpy >0.5, and scipy >0.5. Forward modeling software is also a prerequiste; see _Solver Configuration_ for details.


Hardware Prerequisties
----------------------

Access to a computer cluster is required for most jobs.  Base classes are provided for several common cluster configurations, including PBS and SLURM.  Nonstandard configurations can be accommodated by modifying base classes; see _System Configuration_ for more information.



SOLVER CONFIGURATION
====================

SeisFlows includes python interfaces for SPECFEM2D, SPECFEM3D, and SPECFEM3D_GLOBE.  While these interfaces are part of the main package, the solver source code itself is located elsewhere and must be downloaded separately.  After downloading, users must configure and compile the source code, following the detailed instructions included in the solver user manual.

Summarized briefly, the configuration and compilation procedure is as follows.

Prior to compilation, users need to run the configure script and prepare input files such as

- parameter file
- source file
- stations file.

The result of compilation is a set of binaries files, including

- mesher binary
- solver binary.

Note that if input files change, binary files may need to be recompiled.


Solver Integration
------------------

Integration of the solver with the other workflow components can be challenging. Here, we try to give an idea of the issues involved from both a developer and a user standpoint.

First, developers must decide how to approach compilation.  Because solver computations account for a large portion of the total cost of an inversion, the solver is usually written in an efficient compiled language. If the compilation procedure is not too complex, it may be possible to provde some mechanism for automatically compiling solver binaries. Unfortunately, because the compilation procedure often changes from version to version, there is currently no such mechanism for compiling SPECFEM2D, SPECFEM3D, or SPECFEM3D_GLOBE source code. Users must supply there own binaries for these solvers.

A second issue involves argument passing.

A third issue involves data types. In the solver interface class, it is natural to store models as dictionaries, with different keys corresponding to different material parameters.


Writing Custom Solver Interfaces
--------------------------------

Besides SPECFEM2D, SPECFEM3D, and SPECFEM3D_GLOBE, SeisFlows can interface with other solvers capable of running forward and adjoint simulations. For information about writing custom solver interfaces, see _Developer Reference_.



SYSTEM CONFIGURATION
====================

SeisFlows can run on SLURM, PBS Torque, and PBS PRO clusters.

Job submission on PBS and SLURM systems shares many of the same generalities, but some of the specifics are different.  Our approach to such differences is to create a thin python layer over system commands. For example, system specific job submission commands, such as qsub on PBS or sbatch on SLURM, are replace with a consistent python interface.

If a forward simulation spans one node or less, users should select 'pbs' or 'slurm'. If a forward simulation spans more than one node, users should select 'pbs_big_job' or 'slurm_big_job'.

For debugging, an option to run simulations in serial is also provided, ...



JOB SUBMISSION
==============

Each job must be submitted from its own working directory.  Within a working directory, users must supply two input files, paths.py and paramters.py, which are described below. Output files, by default, are written to the working directory, along with scratch files created by the solver and optimization routines. Different output and scratch directories can be specified by adding or modifying entries in paths.py.

parameters.py
Contains a list of parameter names and values. Prior to a job being submitted, parameters are checked so that errors can be detected without loss of queuetime or walltime. Parameters are stored in dictionary that is accessible from anywhere in the python code. By convention, all parameter names must be upper case. Parameter values can be floats, integers, strings or any other python data type. Parameters can be listed in any order.

paths.py
Contains a list of path names and values. Prior to a job being submitted, paths are checked so that errors can be detected without loss of queuetime or walltime. Paths are stored in a dictionary that is is accessible from anywhere in the python code. By convention, all names must be upper case, and all values must be absolute paths. Paths can be listed in any order.

Once a working directory and input files have been created, users can type sfrun from within the working directory to submit a job. If the 'serial' system configuration is specified in parameters.py, the job will begin executing immediately. If 'pbs' or 'slurm' system configurations are specified, the job will run when resources become available. Once the job begins running, status information will be displayed either to the terminal or to the file ``output.log``.



DEVELOPER REFERENCE
===================

To allow classes to work with one another, each class must conform to an established interface.  In practice, this means that each class must implement a specificed set of methods, with specified input and output.

Methods expected for each class are listed below. Besides required methods, classes may include any number of private methods or utility functions.

