2D Example Walkthrough
======================

The notebook below details a walkthrough of Example \#1 shown in the `SPECFEM2D example \#1 <specfem2d_example.html>`__. This is meant for those who want to understand what is going on under the hood. You are welcome to follow along on your workstation. The following Table of Contents outlines the steps we will take in this tutorial:

.. warning:: 
    Navigation links will not work outside of Jupyter. Please use the navigation bar to the left.

1. `Setup SPECFEM2D <#1.-Setup-SPECFEM2D>`__

   a. `Download and compile
      codebase <#1a.-Download-and-compile-codebase*>`__
   b. `Create a separate SPECFEM2D working
      directory <#1b.-Create-a-separate-SPECFEM2D-working-directory>`__
   c. `Generate initial and target
      models <#1c.-Generate-initial-and-target-models>`__

2. `Initialize SeisFlows (SF) <#2.-Initialize-SeisFlows-(SF)>`__

   a. `SeisFlows working directory and parameter
      file <#2a.-SF-working-directory-and-parameter-file>`__

3. `Run SeisFlows <#2.-Run-SeisFlows>`__

   a. `Forward simulations <#3a.-Forward-simulations>`__
   b. `Exploring the SeisFlows directory
      structure <#3b.-Exploring-the-SF-directory-structure>`__
   c. `Adjoint simulations <#3c.-Adjoint-simulations>`__
   d. `Line search and model
      update <#3d.-Line-search-and-model-update>`__

4. `Conclusions <#4.-Conclusions>`__

1. Setup SPECFEM2D
~~~~~~~~~~~~~~~~~~

1a. Download and compile codebase (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   **NOTE**: If you have already downloaded and compiled SPECFEM2D, you
   can skip most of this subsection (1a). However you will need to edit
   the first two paths in the following cell (WORKDIR and
   SPECFEM2D_ORIGINAL), and execute the path structure defined in the
   cell.

First we’ll download and compile SPECFEM2D to generate the binaries
necessary to run our simulations. We will then populate a new SPECFEM2D
working directory that will be used by SeisFlows. We’ll use to Python OS
module to do our filesystem processes just to keep everything in Python,
but this can easily be accomplished in bash.

.. code:: ipython3

    import os
    import glob
    import shutil
    import numpy as np

.. code:: ipython3

    # vvv USER MUST EDIT THE FOLLOWING PATHS vvv
    WORKDIR = "/home/bchow/Work/scratch" 
    SPECFEM2D = "/home/bchow/REPOSITORIES/specfem2d"
    # where WORKDIR: points to your own working directory
    # and SPECFEM2D: points to an existing specfem2D repository if available (if not set as '')
    # ^^^ USER MUST EDIT THE FOLLOWING PATHS ^^^
    # ======================================================================================================
    
    # Distribute the necessary file structure of the SPECFEM2D repository that we will downloaded/reference
    SPECFEM2D_ORIGINAL = os.path.join(WORKDIR, "specfem2d")
    SPECFEM2D_BIN_ORIGINAL = os.path.join(SPECFEM2D_ORIGINAL, "bin")
    SPECFEM2D_DATA_ORIGINAL = os.path.join(SPECFEM2D_ORIGINAL, "DATA")
    TAPE_2007_EXAMPLE = os.path.join(SPECFEM2D_ORIGINAL, "EXAMPLES", "Tape2007")
    
    # The SPECFEM2D working directory that we will create separate from the downloaded repo
    SPECFEM2D_WORKDIR = os.path.join(WORKDIR, "specfem2d_workdir")
    SPECFEM2D_BIN = os.path.join(SPECFEM2D_WORKDIR, "bin")
    SPECFEM2D_DATA = os.path.join(SPECFEM2D_WORKDIR, "DATA")
    SPECFEM2D_OUTPUT = os.path.join(SPECFEM2D_WORKDIR, "OUTPUT_FILES")
    
    # Pre-defined locations of velocity models we will generate using the solver
    SPECFEM2D_MODEL_INIT = os.path.join(SPECFEM2D_WORKDIR, "OUTPUT_FILES_INIT")
    SPECFEM2D_MODEL_TRUE = os.path.join(SPECFEM2D_WORKDIR, "OUTPUT_FILES_TRUE")

.. code:: ipython3

    # Download SPECFEM2D from GitHub, devel branch for latest codebase OR symlink from existing repo
    if not os.path.exists(WORKDIR):
        os.makedirs(WORKDIR)
    os.chdir(WORKDIR)
    
    if os.path.exists("specfem2d"):
        print("SPECFEM2D repository already found, you may skip this subsection")
        pass
    elif os.path.exists(SPECFEM2D):
        print("Existing SPECMFE2D respository found, symlinking to working directory")
        os.symlink(SPECFEM2D, "./specfem2d")
    else:
        print("Cloning respository from GitHub")
        ! git clone --recursive --branch devel https://github.com/geodynamics/specfem2d.git


.. parsed-literal::

    Existing SPECMFE2D respository found, symlinking to working directory


.. code:: ipython3

    # Compile SPECFEM2D to generate the Makefile
    os.chdir(SPECFEM2D_ORIGINAL)
    if not os.path.exists("./config.log"):
        os.system("./configure")

.. code:: ipython3

    # Run make to generate SPECFEM2D binaries
    if not os.path.exists("bin"):
        os.system("make all")

.. code:: ipython3

    # Check out the binary files that have been created
    os.chdir(SPECFEM2D_ORIGINAL)
    ! pwd
    ! ls bin/


.. parsed-literal::

    /home/bchow/REPOSITORIES/specfem2d
    xadj_seismogram		      xconvolve_source_timefunction  xspecfem2D
    xcheck_quality_external_mesh  xmeshfem2D		     xsum_kernels
    xcombine_sem		      xsmooth_sem


1b. Create a separate SPECFEM2D working directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next we’ll create a new SPECFEM2D working directory, separate from the
original repository. The intent here is to isolate the original
SPECFEM2D repository from our working state, to protect it from things
like accidental file deletions or manipulations. This is not a mandatory
step for using SeisFlows, but it helps keep file structure clean in the
long run, and is the SeisFlows3 dev team’s preferred method of using
SPECFEM.

.. note::
    All SPECFEM2D/3D/3D_GLOBE need to run successfully are the bin/, DATA/, and OUTPUT_FILES/ directories. Everything else in the repository is not mandatory for running binaries.

In this tutorial we will be using the `Tape2007 example
problem <https://github.com/geodynamics/specfem2d/tree/devel/EXAMPLES/Tape2007>`__
to define our **DATA/** directory (last tested 8/15/22, bdba4389).

.. code:: ipython3

    # Incase we've run this docs page before, delete the working directory before remaking
    if os.path.exists(SPECFEM2D_WORKDIR):
        shutil.rmtree(SPECFEM2D_WORKDIR)
    
    os.mkdir(SPECFEM2D_WORKDIR)
    os.chdir(SPECFEM2D_WORKDIR)
    
    # Copy the binary files incase we update the source code. These can also be symlinked.
    shutil.copytree(SPECFEM2D_BIN_ORIGINAL, "bin")
    
    # Copy the DATA/ directory because we will be making edits here frequently and it's useful to
    # retain the original files for reference. We will be running one of the example problems: Tape2007
    shutil.copytree(os.path.join(TAPE_2007_EXAMPLE, "DATA"), "DATA")
    
    ! pwd
    ! ls


.. parsed-literal::

    /home/bchow/Work/scratch/specfem2d_workdir
    bin  DATA


.. code:: ipython3

    # Run the Tape2007 example to make sure SPECFEM2D is working as expected
    os.chdir(TAPE_2007_EXAMPLE)
    ! ./run_this_example.sh > output_log.txt
    
    assert(os.path.exists("OUTPUT_FILES/forward_image000004800.jpg")), \
        (f"Example did not run, the remainder of this docs page will likely not work."
         f"Please check the following directory: {TAPE_2007_EXAMPLE}")
    
    ! tail output_log.txt


.. parsed-literal::

     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------
     D a t e : 16 - 08 - 2022                                 T i m e  : 14:26:37
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------
    
    see results in directory: OUTPUT_FILES/
    
    done
    Tue Aug 16 02:26:37 PM AKDT 2022


--------------

Now we need to manually set up our SPECFEM2D working directory. As
mentioned in the previous cell, the only required elements of this
working directory are the following (these files will form the basis for
how SeisFlows3 operates within the SPECFEM2D framework):

1. **bin/** directory containing SPECFEM2D binaries
2. **DATA/** directory containing SOURCE and STATION files, as well as a
   SPECFEM2D Par_file
3. \__OUTPUT_FILES/proc??????_*.bin_\_ files which define the starting
   (and target) models

.. note:: 
    This file structure is the same for all versions of SPECFEM (2D/3D/3D_GLOBE)

.. code:: ipython3

    # First we will set the correct SOURCE and STATION files.
    # This is the same task as shown in ./run_this_example.sh
    os.chdir(SPECFEM2D_DATA)
    
    # Symlink source 001 as our main source
    if os.path.exists("SOURCE"):
        os.remove("SOURCE")
    os.symlink("SOURCE_001", "SOURCE")
    
    # Copy the correct Par_file so that edits do not affect the original file
    if os.path.exists("Par_file"):
        os.remove("Par_file")
    shutil.copy("Par_file_Tape2007_onerec", "Par_file")
    
    ! ls


.. parsed-literal::

    interfaces_Tape2007.dat		     SOURCE_003  SOURCE_012  SOURCE_021
    model_velocity.dat_checker	     SOURCE_004  SOURCE_013  SOURCE_022
    Par_file			     SOURCE_005  SOURCE_014  SOURCE_023
    Par_file_Tape2007_132rec_checker     SOURCE_006  SOURCE_015  SOURCE_024
    Par_file_Tape2007_onerec	     SOURCE_007  SOURCE_016  SOURCE_025
    proc000000_model_velocity.dat_input  SOURCE_008  SOURCE_017  STATIONS
    SOURCE				     SOURCE_009  SOURCE_018  STATIONS_checker
    SOURCE_001			     SOURCE_010  SOURCE_019
    SOURCE_002			     SOURCE_011  SOURCE_020


1c. Generate initial and target models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since we’re doing a synthetic-synthetic inversion, we need to manually
set up the velocity models with which we generate our synthetic
waveforms. The naming conventions for these models are:

1. **MODEL_INIT:** The initial or starting model. Used to generate the
   actual synthetic seismograms. This is considered M00.
2. **MODEL_TRUE:** The target or true model. Used to generate ‘data’
   (also synthetic). This is the reference model that our inversion is
   trying to resolve.

The starting model is defined as a homogeneous halfspace uin the
Tape2007 example problem. We will need to run both ``xmeshfem2D`` and
``xspecfem2D`` to generate the required velocity model database files.
We will generate our target model by slightly perturbing the parameters
of the initial model.

.. note::
    We can use the SeisFlows3 command line option `seisflows sempar` to directly edit the SPECFEM2D Par_file in the command line. This will work for the SPECFEM3D Par_file as well.

.. code:: ipython3

    os.chdir(SPECFEM2D_DATA)
    
    # Ensure that SPECFEM2D outputs the velocity model in the expected binary format
    ! seisflows sempar setup_with_binary_database 1  # allow creation of .bin files
    ! seisflows sempar save_model binary  # output model in .bin database format
    ! seisflows sempar save_ascii_kernels .false.  # output kernels in .bin format, not ASCII


.. parsed-literal::

    setup_with_binary_database: 0 -> 1
    SAVE_MODEL: default -> binary
    save_ASCII_kernels: .true. -> .false.


.. code:: ipython3

    # SPECFEM requires that we create the OUTPUT_FILES directory before running
    os.chdir(SPECFEM2D_WORKDIR)
    
    if os.path.exists(SPECFEM2D_OUTPUT):
        shutil.rmtree(SPECFEM2D_OUTPUT)
        
    os.mkdir(SPECFEM2D_OUTPUT)
    
    ! ls


.. parsed-literal::

    bin  DATA  OUTPUT_FILES


.. code:: ipython3

    # GENERATE MODEL_INIT
    os.chdir(SPECFEM2D_WORKDIR)
    
    # Run the mesher and solver to generate our initial model
    ! ./bin/xmeshfem2D > OUTPUT_FILES/mesher_log.txt
    ! ./bin/xspecfem2D > OUTPUT_FILES/solver_log.txt
    
    # Move the model files (*.bin) into the OUTPUT_FILES directory, where SeisFlows3 expects them
    ! mv DATA/*bin OUTPUT_FILES
    
    # Make sure we don't overwrite this initial model when creating our target model in the next step
    ! mv OUTPUT_FILES OUTPUT_FILES_INIT
    
    ! head OUTPUT_FILES_INIT/solver_log.txt
    ! tail OUTPUT_FILES_INIT/solver_log.txt


.. parsed-literal::

    
     **********************************************
     **** Specfem 2-D Solver - serial version  ****
     **********************************************
    
     Running Git version of the code corresponding to commit cf89366717d9435985ba852ef1d41a10cee97884
     dating From Date:   Mon Nov 29 23:20:51 2021 -0800
    
    
     NDIM =            2
     -------------------------------------------------------------------------------
     Program SPECFEM2D: 
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------
     Tape-Liu-Tromp (GJI 2007)
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------
     D a t e : 16 - 08 - 2022                                 T i m e  : 14:26:52
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------


--------------

Now we want to perturb the initial model to create our target model
(**MODEL_TRUE**). The seisflows command line subargument
``seisflows sempar velocity_model`` will let us view and edit the
velocity model. You can also do this manually by editing the Par_file
directly.

.. code:: ipython3

    # GENERATE MODEL_TRUE
    os.chdir(SPECFEM2D_DATA)
    
    # Edit the Par_file by increasing velocities by ~10% 
    ! seisflows sempar velocity_model '1 1 2600.d0 5900.d0 3550.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0'


.. parsed-literal::

    VELOCITY_MODEL:
    
    1 1 2600.d0 5800.d0 3500.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0
    ->
    1 1 2600.d0 5900.d0 3550.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0


.. code:: ipython3

    # Re-run the mesher and solver to generate our target velocity model
    os.chdir(SPECFEM2D_WORKDIR)
    
    # Make sure the ./OUTPUT_FILES directory exists since we moved the old one
    if os.path.exists(SPECFEM2D_OUTPUT):
        shutil.rmtree(SPECFEM2D_OUTPUT)
    os.mkdir(SPECFEM2D_OUTPUT)
    
    # Run the binaries to generate MODEL_TRUE
    ! ./bin/xmeshfem2D > OUTPUT_FILES/mesher_log.txt
    ! ./bin/xspecfem2D > OUTPUT_FILES/solver_log.txt
    
    # Move all the relevant files into OUTPUT_FILES 
    ! mv ./DATA/*bin OUTPUT_FILES
    ! mv OUTPUT_FILES OUTPUT_FILES_TRUE
    
    ! head OUTPUT_FILES_INIT/solver_log.txt
    ! tail OUTPUT_FILES_INIT/solver_log.txt


.. parsed-literal::

    
     **********************************************
     **** Specfem 2-D Solver - serial version  ****
     **********************************************
    
     Running Git version of the code corresponding to commit cf89366717d9435985ba852ef1d41a10cee97884
     dating From Date:   Mon Nov 29 23:20:51 2021 -0800
    
    
     NDIM =            2
     -------------------------------------------------------------------------------
     Program SPECFEM2D: 
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------
     Tape-Liu-Tromp (GJI 2007)
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------
     D a t e : 16 - 08 - 2022                                 T i m e  : 14:26:52
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------


.. code:: ipython3

    # Great, we have all the necessary SPECFEM files to run our SeisFlows inversion!
    ! ls


.. parsed-literal::

    bin  DATA  OUTPUT_FILES_INIT  OUTPUT_FILES_TRUE


2. Initialize SeisFlows (SF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this Section we will look at a SeisFlows working directory, parameter
file, and working state.

2a. SeisFlows working directory and parameter file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As with SPECFEM, SeisFlows requires a parameter file
(**parameters.yaml**) that controls how an automated workflow will
proceed. Because SeisFlows is modular, there are a large number of
potential parameters which may be present in a SeisFlows parameter file,
as each sub-module may have its own set of unique parameters.

In contrast to SPECFEM’s method of listing all available parameters and
leaving it up the User to determine which ones are relevant to them,
SeisFlows dynamically builds its parameter file based on User inputs. In
this subsection we will use the built-in SeisFlows command line tools to
generate and populate the parameter file.

.. note::
    See the `parameter file documentation page <parameter_file.html>`__ for a more in depth exploration of this central SeisFlows file.

In the previous section we saw the ``sempar`` command in action. We can
use the ``-h`` or help flag to list all available SiesFlows3 command
line commands.

.. code:: ipython3

    ! seisflows -h


.. parsed-literal::

    usage: seisflows [-h] [-w [WORKDIR]] [-p [PARAMETER_FILE]]
                     {setup,configure,swap,init,submit,resume,restart,clean,par,sempar,check,print,reset,debug,examples}
                     ...
    
    ================================================================================
    
                         SeisFlows: Waveform Inversion Package                      
    
    ================================================================================
    
    optional arguments:
      -h, --help            show this help message and exit
      -w [WORKDIR], --workdir [WORKDIR]
                            The SeisFlows working directory, default: cwd
      -p [PARAMETER_FILE], --parameter_file [PARAMETER_FILE]
                            Parameters file, default: 'parameters.yaml'
    
    command:
      Available SeisFlows arguments and their intended usages
    
        setup               Setup working directory from scratch
        configure           Fill parameter file with defaults
        swap                Swap module parameters in an existing parameter file
        init                Initiate working environment
        submit              Submit initial workflow to system
        resume              Re-submit previous workflow to system
        restart             Remove current environment and submit new workflow
        clean               Remove files relating to an active working environment
        par                 View and edit SeisFlows parameter file
        sempar              View and edit SPECFEM parameter file
        check               Check state of an active environment
        print               Print information related to an active environment
        reset               Reset modules within an active state
        debug               Start interactive debug environment
        examples            Look at and run pre-configured example problems
    
    'seisflows [command] -h' for more detailed descriptions of each command.


.. code:: ipython3

    # The command 'setup' creates the 'parameters.yaml' file that controls all of SeisFlows
    # the '-f' flag removes any exist 'parameters.yaml' file that might be in the directory
    os.chdir(WORKDIR)
    ! seisflows setup -f
    ! ls


.. parsed-literal::

    creating parameter file: parameters.yaml
    parameters.yaml  sflog.txt  specfem2d  specfem2d_workdir


.. code:: ipython3

    # Let's have a look at this file, which has not yet been populated
    ! cat parameters.yaml


.. parsed-literal::

    # //////////////////////////////////////////////////////////////////////////////
    #
    #                        SeisFlows YAML Parameter File
    #
    # //////////////////////////////////////////////////////////////////////////////
    #
    # Modules correspond to the structure of the source code, and determine
    # SeisFlows' behavior at runtime. Each module requires its own sub-parameters.
    #
    # .. rubric::
    #   - To determine available options for modules listed below, run:
    #       > seisflows print modules
    #   - To auto-fill with docstrings and default values (recommended), run:
    #       > seisflows configure
    #   - To set values as NoneType, use: null
    #   - To set values as infinity, use: inf
    #
    #                                    MODULES
    #                                    ///////
    # workflow (str):    The types and order of functions for running SeisFlows
    # system (str):      Computer architecture of the system being used
    # solver (str):      External numerical solver to use for waveform simulations
    # preprocess (str):  Preprocessing schema for waveform data
    # optimize (str):    Optimization algorithm for the inverse problem
    # ==============================================================================
    workflow: forward
    system: workstation
    solver: specfem2d
    preprocess: default
    optimize: gradient


.. code:: ipython3

    # We can use the `seisflows print modules` command to list out the available options 
    ! seisflows print modules


.. parsed-literal::

                                   SEISFLOWS MODULES                                
                                   /////////////////                                
    '-': module, '*': class
    
    - workflow
        * forward
        * inversion
        * migration
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

    # For this example, we can use most of the default modules, however we need to 
    # change the SOLVER module to let SeisFlows know we're using SPECFEM2D (as opposed to 3D)
    ! seisflows par workflow inversion
    ! cat parameters.yaml


.. parsed-literal::

    workflow: forward -> inversion
    # //////////////////////////////////////////////////////////////////////////////
    #
    #                        SeisFlows YAML Parameter File
    #
    # //////////////////////////////////////////////////////////////////////////////
    #
    # Modules correspond to the structure of the source code, and determine
    # SeisFlows' behavior at runtime. Each module requires its own sub-parameters.
    #
    # .. rubric::
    #   - To determine available options for modules listed below, run:
    #       > seisflows print modules
    #   - To auto-fill with docstrings and default values (recommended), run:
    #       > seisflows configure
    #   - To set values as NoneType, use: null
    #   - To set values as infinity, use: inf
    #
    #                                    MODULES
    #                                    ///////
    # workflow (str):    The types and order of functions for running SeisFlows
    # system (str):      Computer architecture of the system being used
    # solver (str):      External numerical solver to use for waveform simulations
    # preprocess (str):  Preprocessing schema for waveform data
    # optimize (str):    Optimization algorithm for the inverse problem
    # ==============================================================================
    workflow: inversion
    system: workstation
    solver: specfem2d
    preprocess: default
    optimize: gradient


--------------

The ``seisflows configure`` command populates the parameter file based
on the chosen modules. SeisFlows will attempt to fill in all parameters
with reasonable default values. Docstrings above each module show
descriptions and available options for each of these parameters.

In the follownig cell we will use the ``seisflows par`` command to edit
the parameters.yaml file directly, replacing some default parameters
with our own values. Comments next to each evaluation describe the
choice for each.

.. code:: ipython3

    ! seisflows configure
    ! head --lines=50 parameters.yaml


.. parsed-literal::

    # //////////////////////////////////////////////////////////////////////////////
    #
    #                        SeisFlows YAML Parameter File
    #
    # //////////////////////////////////////////////////////////////////////////////
    #
    # Modules correspond to the structure of the source code, and determine
    # SeisFlows' behavior at runtime. Each module requires its own sub-parameters.
    #
    # .. rubric::
    #   - To determine available options for modules listed below, run:
    #       > seisflows print modules
    #   - To auto-fill with docstrings and default values (recommended), run:
    #       > seisflows configure
    #   - To set values as NoneType, use: null
    #   - To set values as infinity, use: inf
    #
    #                                    MODULES
    #                                    ///////
    # workflow (str):    The types and order of functions for running SeisFlows
    # system (str):      Computer architecture of the system being used
    # solver (str):      External numerical solver to use for waveform simulations
    # preprocess (str):  Preprocessing schema for waveform data
    # optimize (str):    Optimization algorithm for the inverse problem
    # ==============================================================================
    workflow: inversion
    system: workstation
    solver: specfem2d
    preprocess: default
    optimize: gradient
    # =============================================================================
    #
    #    Forward Workflow
    #    ----------------
    #    Run forward solver in parallel and (optionally) calculate
    #    data-synthetic misfit and adjoint sources.
    #
    #    Parameters
    #    ----------
    #    :type modules: list of module
    #    :param modules: instantiated SeisFlows modules which should have been
    #        generated by the function `seisflows.config.import_seisflows` with a
    #        parameter file generated by seisflows.configure
    #    :type data_case: str
    #    :param data_case: How to address 'data' in the workflow, available options:
    #        'data': real data will be provided by the user in
    #        `path_data/{source_name}` in the same format that the solver will
    #        produce synthetics (controlled by `solver.format`) OR
    #        synthetic': 'data' will be generated as synthetic seismograms using
    #        a target model provided in `path_model_true`. If None, workflow will


.. code:: ipython3

    # EDIT THE SEISFLOWS PARAMETER FILE
    ! seisflows par ntask 3  # set the number of sources/events to use
    ! seisflows par materials elastic  # update Vp and Vs during inversion
    ! seisflows par end 2  # final iteration -- we will only run 1
    ! seisflows par data_case synthetic  # synthetic-synthetic means we need both INIT and TRUE models
    ! seisflows par components Y  # this default example creates Y-component seismograms
    ! seisflows par step_count_max 5  # limit the number of steps in the line search
    
    # Use Python syntax here to access path constants
    os.system(f"seisflows par path_specfem_bin {SPECFEM2D_BIN}")  # set path to SPECFEM2D binaries
    os.system(f"seisflows par path_specfem_data {SPECFEM2D_DATA}")  # set path to SEPCFEM2D DATA/
    os.system(f"seisflows par path_model_init {SPECFEM2D_MODEL_INIT}")  # set path to INIT model
    os.system(f"seisflows par path_model_true {SPECFEM2D_MODEL_TRUE}")  # set path to TRUE model


.. parsed-literal::

    ntask: 1 -> 3
    materials: acoustic -> elastic
    end: 1 -> 2
    data_case: data -> synthetic
    components: ZNE -> Y
    step_count_max: 10 -> 5
    path_specfem_bin: null -> /home/bchow/Work/scratch/specfem2d_workdir/bin
    path_specfem_data: null -> /home/bchow/Work/scratch/specfem2d_workdir/DATA
    path_model_init: null ->
    /home/bchow/Work/scratch/specfem2d_workdir/OUTPUT_FILES_INIT
    path_model_true: null ->
    /home/bchow/Work/scratch/specfem2d_workdir/OUTPUT_FILES_TRUE




.. parsed-literal::

    0



--------------

One last thing, we will need to edit the SPECFEM2D Par_file parameter
``MODEL`` such that ``xmeshfem2d`` reads our pre-built velocity models
(*.bin files) rather than the meshing parameters defined in the
Par_file.

.. code:: ipython3

    os.chdir(SPECFEM2D_DATA)
    ! seisflows sempar model gll


.. parsed-literal::

    MODEL: default -> gll


3. Run SeisFlows
~~~~~~~~~~~~~~~~

In this Section we will run SeisFlows to generate synthetic seismograms,
kernels, a gradient, and an updated velocity model.

3a. Forward simulations
^^^^^^^^^^^^^^^^^^^^^^^

SeisFlows is an automated workflow tool, such that once we run
``seisflows submit`` we should not need to intervene in the workflow.
However the package does allow the User flexibility in how they want the
workflow to behave.

For example, we can run our workflow in stages by taking advantage of
the ``stop_after`` parameter. As its name suggests, ``stop_after``
allows us to stop a workflow prematurely so that we may stop and look at
results, or debug a failing workflow.

The ``seisflows print flow`` command tells us what functions we can use
for the ``stop_after`` parameter.

.. code:: ipython3

    os.chdir(WORKDIR)
    ! seisflows print tasks


.. parsed-literal::

                              SEISFLOWS WORKFLOW TASK LIST                          
                              ////////////////////////////                          
    Task list for <class 'seisflows.workflow.inversion.Inversion'>
    
    1: evaluate_initial_misfit
    2: run_adjoint_simulations
    3: postprocess_event_kernels
    4: evaluate_gradient_from_kernels
    5: initialize_line_search
    6: perform_line_search
    7: finalize_iteration


--------------

In the Inversion workflow, the tasks listed are described as follows:

1. **evaluate_initial_misfit:**

   a. Prepare data for inversion by either copying data from disk or
      generating ‘synthetic data’ with MODEL_TRUE
   b. Call numerical solver to run forward simulations using MODEL_INIT,
      generating synthetics
   c. Evaluate the objective function by performing waveform comparisons
   d. Prepare ``run_adjoint_simulations`` step by generating adjoint
      sources and auxiliary files

2. **run_adjoint_simulations:** Call numerical solver to run adjoint
   simulation, generating kernels
3. **postprocess_event_kernels:** Combine all event kernels into a
   misfit kernel.
4. **evaluate_gradient_from_kernels:** Smooth and mask the misfit kernel
   to create the gradient
5. **initialize_line_search:** Call on the optimization library to scale
   the gradient by a step length to compute the search direction.
   Prepare file structure for line search.
6. **perform_line_search:** Perform a line search by algorithmically
   scaling the gradient and evaluating the misfit function (forward
   simulations and misfit quantification) until misfit is acceptably
   reduced.
7. **finalize_iteration:** Run any finalization steps such as saving
   traces, kernels, gradients and models to disk, setting up SeisFlows3
   for any subsequent iterations. Clean the scratch/ directory in
   preparation for subsequent iterations

Let’s set the ``stop_after`` argument to **evaluate_initial_misfit**,
this will halt the workflow after the intialization step.

.. code:: ipython3

    ! seisflows par stop_after evaluate_initial_misfit


.. parsed-literal::

    stop_after: null -> evaluate_initial_misfit


--------------

Now let’s run SeisFlows. There are two ways to do this: ``submit`` and
``restart``

1. ``seisflows submit`` is used to run new workflows and resume stopped
   or failed workflows.
2. The ``restart`` command is simply a convenience function that runs
   ``clean`` (to remove an active working state) and ``submit`` (to
   submit a fresh workflow).

Since this is our first run, we’ll use ``seisflows submit``.

.. code:: ipython3

    ! seisflows submit 


.. parsed-literal::

    2022-08-16 14:32:48 (I) | 
    ================================================================================
                             SETTING UP INVERSION WORKFLOW                          
    ================================================================================
    2022-08-16 14:32:55 (D) | running setup for module 'system.Workstation'
    2022-08-16 14:32:57 (D) | copying par/log file to: /home/bchow/Work/scratch/logs/sflog_001.txt
    2022-08-16 14:32:57 (D) | copying par/log file to: /home/bchow/Work/scratch/logs/parameters_001.yaml
    2022-08-16 14:32:57 (D) | running setup for module 'solver.Specfem2D'
    2022-08-16 14:32:57 (I) | initializing 3 solver directories
    2022-08-16 14:32:57 (D) | initializing solver directory source: 001
    2022-08-16 14:33:04 (D) | linking source '001' as 'mainsolver'
    2022-08-16 14:33:04 (D) | initializing solver directory source: 002
    2022-08-16 14:33:09 (D) | initializing solver directory source: 003
    2022-08-16 14:33:16 (D) | running setup for module 'preprocess.Default'
    2022-08-16 14:33:16 (D) | running setup for module 'optimize.Gradient'
    2022-08-16 14:33:17 (I) | no optimization checkpoint found, assuming first run
    2022-08-16 14:33:17 (I) | re-loading optimization module from checkpoint
    2022-08-16 14:33:17 (I) | 
    ////////////////////////////////////////////////////////////////////////////////
                                  RUNNING ITERATION 01                              
    ////////////////////////////////////////////////////////////////////////////////
    2022-08-16 14:33:17 (I) | 
    ================================================================================
                               RUNNING INVERSION WORKFLOW                           
    ================================================================================
    2022-08-16 14:33:17 (I) | 
    ////////////////////////////////////////////////////////////////////////////////
                          EVALUATING MISFIT FOR INITIAL MODEL                       
    ////////////////////////////////////////////////////////////////////////////////
    2022-08-16 14:33:17 (I) | checking initial model parameters
    2022-08-16 14:33:17 (I) | 5800.00 <= vp <= 5800.00
    2022-08-16 14:33:17 (I) | 2600.00 <= rho <= 2600.00
    2022-08-16 14:33:17 (I) | 3500.00 <= vs <= 3500.00
    2022-08-16 14:33:17 (I) | checking true/target model parameters
    2022-08-16 14:33:17 (I) | 5900.00 <= vp <= 5900.00
    2022-08-16 14:33:17 (I) | 2600.00 <= rho <= 2600.00
    2022-08-16 14:33:17 (I) | 3550.00 <= vs <= 3550.00
    2022-08-16 14:33:17 (I) | preparing observation data for source 001
    2022-08-16 14:33:17 (I) | running forward simulation w/ target model for 001
    2022-08-16 14:33:21 (I) | evaluating objective function for source 001
    2022-08-16 14:33:21 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:33:25 (D) | quantifying misfit with 'Default'
    2022-08-16 14:33:25 (I) | preparing observation data for source 002
    2022-08-16 14:33:25 (I) | running forward simulation w/ target model for 002
    2022-08-16 14:33:29 (I) | evaluating objective function for source 002
    2022-08-16 14:33:29 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:33:33 (D) | quantifying misfit with 'Default'
    2022-08-16 14:33:33 (I) | preparing observation data for source 003
    2022-08-16 14:33:33 (I) | running forward simulation w/ target model for 003
    2022-08-16 14:33:36 (I) | evaluating objective function for source 003
    2022-08-16 14:33:36 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:33:40 (D) | quantifying misfit with 'Default'
    2022-08-16 14:33:40 (I) | stop workflow at `stop_after`: evaluate_initial_misfit


.. note::
    For a detailed exploration of a SeisFlows working directory, see the `working directory <working_directory.html>`__ documentation page where we explain each of the files and directories that have been generated during this workflow. Below we just look at two files which are required for our adjoint simulation, the adjoint sources (.adj) and STATIONS_ADJOINT file

.. code:: ipython3

    # The adjoint source is created in the same format as the synthetics (two-column ASCII) 
    ! head scratch/solver/001/traces/adj/AA.S0001.BXY.adj


.. parsed-literal::

      -48.0000000         0.0000000
      -47.9400000         0.0000000
      -47.8800000         0.0000000
      -47.8200000         0.0000000
      -47.7600000         0.0000000
      -47.7000000         0.0000000
      -47.6400000         0.0000000
      -47.5800000         0.0000000
      -47.5200000         0.0000000
      -47.4600000         0.0000000


3b. Adjoint simulations
^^^^^^^^^^^^^^^^^^^^^^^

Now that we have all the required files for running an adjoint
simulation (*.adj waveforms and STATIONS_ADJOINT file), we can continue
with the SeisFlows3 Inversion workflow. No need to edit the Par_file or
anything like that, SeisFlows3 will take care of that under the hood. We
simply need to tell the workflow (via the parameters.yaml file) to
``resume_from`` the correct function. We can have a look at these
functions again:

.. code:: ipython3

    ! seisflows print tasks


.. parsed-literal::

                              SEISFLOWS WORKFLOW TASK LIST                          
                              ////////////////////////////                          
    Task list for <class 'seisflows.workflow.inversion.Inversion'>
    
    1: evaluate_initial_misfit
    2: run_adjoint_simulations
    3: postprocess_event_kernels
    4: evaluate_gradient_from_kernels
    5: initialize_line_search
    6: perform_line_search
    7: finalize_iteration


.. code:: ipython3

    # We'll stop just before the line search so that we can take a look at the files 
    # generated during the middle tasks
    ! seisflows par stop_after evaluate_gradient_from_kernels


.. parsed-literal::

    stop_after: evaluate_initial_misfit -> evaluate_gradient_from_kernels


.. code:: ipython3

    # We can use the `seisflows submit` command to continue an active workflow
    # The state file created during the first run will tell the workflow to resume from the stopped point in the workflow
    ! seisflows submit 


.. parsed-literal::

    2022-08-16 14:36:42 (D) | setting iteration==1 from state file
    2022-08-16 14:36:42 (I) | 
    ================================================================================
                             SETTING UP INVERSION WORKFLOW                          
    ================================================================================
    2022-08-16 14:36:48 (D) | running setup for module 'system.Workstation'
    2022-08-16 14:36:51 (D) | copying par/log file to: /home/bchow/Work/scratch/logs/sflog_002.txt
    2022-08-16 14:36:51 (D) | copying par/log file to: /home/bchow/Work/scratch/logs/parameters_002.yaml
    2022-08-16 14:36:51 (D) | running setup for module 'solver.Specfem2D'
    2022-08-16 14:36:51 (I) | initializing 3 solver directories
    2022-08-16 14:36:51 (D) | running setup for module 'preprocess.Default'
    2022-08-16 14:36:52 (D) | running setup for module 'optimize.Gradient'
    2022-08-16 14:36:53 (I) | re-loading optimization module from checkpoint
    2022-08-16 14:36:54 (I) | re-loading optimization module from checkpoint
    2022-08-16 14:36:54 (I) | 
    ////////////////////////////////////////////////////////////////////////////////
                                  RUNNING ITERATION 01                              
    ////////////////////////////////////////////////////////////////////////////////
    2022-08-16 14:36:54 (I) | 
    ================================================================================
                               RUNNING INVERSION WORKFLOW                           
    ================================================================================
    2022-08-16 14:36:54 (I) | 'evaluate_initial_misfit' has already been run, skipping
    2022-08-16 14:36:54 (I) | 
    ////////////////////////////////////////////////////////////////////////////////
                    EVALUATING EVENT KERNELS W/ ADJOINT SIMULATIONS                 
    ////////////////////////////////////////////////////////////////////////////////
    2022-08-16 14:36:54 (I) | running SPECFEM executable bin/xspecfem2D, log to 'adj_solver.log'
    2022-08-16 14:37:05 (D) | renaming output event kernels: 'alpha' -> 'vp'
    2022-08-16 14:37:05 (D) | renaming output event kernels: 'beta' -> 'vs'
    2022-08-16 14:37:05 (I) | running SPECFEM executable bin/xspecfem2D, log to 'adj_solver.log'
    2022-08-16 14:37:16 (D) | renaming output event kernels: 'alpha' -> 'vp'
    2022-08-16 14:37:16 (D) | renaming output event kernels: 'beta' -> 'vs'
    2022-08-16 14:37:18 (I) | running SPECFEM executable bin/xspecfem2D, log to 'adj_solver.log'
    2022-08-16 14:37:29 (D) | renaming output event kernels: 'alpha' -> 'vp'
    2022-08-16 14:37:29 (D) | renaming output event kernels: 'beta' -> 'vs'
    2022-08-16 14:37:30 (I) | 
    ////////////////////////////////////////////////////////////////////////////////
                          GENERATING/PROCESSING MISFIT KERNEL                       
    ////////////////////////////////////////////////////////////////////////////////
    2022-08-16 14:37:30 (I) | combining event kernels into single misfit kernel
    2022-08-16 14:37:31 (I) | scaling gradient to absolute model perturbations
    2022-08-16 14:37:32 (I) | stop workflow at `stop_after`: evaluate_gradient_from_kernels


--------------

The function **run_adjoint_simulations()** has run adjoint simulations
to generate event kernels. The functions **postprocess_event_kernels**
and **evaluate_gradient_from_kernels** will have summed and (optionally)
smoothed the kernels to recover the gradient, which will be used to
update our starting model.

   **NOTE**: Since we did not specify any smoothing lenghts
   (PAR.SMOOTH_H and PAR.SMOOTH_V), no smoothing of the gradient has
   occurred.

Using the gradient-descent optimization algorithm, SeisFlows will now
compute a search direction that will be used in the line search to
search for a best fitting model which optimally reduces the objective
function. We can take a look at where SeisFlows has stored the
information relating to kernel generation and the optimization
computation.

.. code:: ipython3

    # Gradient evaluation files are stored here, the kernels are stored separately from the gradient incase
    # the user wants to manually manipulate them
    ! ls scratch/eval_grad


.. parsed-literal::

    gradient  kernels  misfit_kernel  model  residuals.txt


.. code:: ipython3

    # SeisFlows3 stores all kernels and gradient information as SPECFEM binary (.bin) files
    ! ls scratch/eval_grad/gradient


.. parsed-literal::

    proc000000_vp_kernel.bin  proc000000_vs_kernel.bin


.. code:: ipython3

    # Kernels are stored on a per-event basis, and summed together (sum/). If smoothing was performed, 
    # we would see both smoothed and unsmoothed versions of the misfit kernel
    ! ls scratch/eval_grad/kernels


.. parsed-literal::

    001  002  003


.. code:: ipython3

    # We can see that some new values have been stored in prepartion for the line search,
    # including g_new (current gradient) and p_new (current search direction). These are also
    # stored as vector NumPy arrays (.npy files)
    ! ls scratch/optimize


.. parsed-literal::

    checkpoint.npz	f_new.txt  g_new.npz  m_new.npz


.. code:: ipython3

    g_new = np.load("scratch/optimize/g_new.npz")
    print(g_new["vs_kernel"])


.. parsed-literal::

    [[-1.18126331e-12  2.40273470e-12  3.97045036e-11 ...  9.62017688e-11
       4.21140102e-11  3.96825021e-12]]


--------------

3c. Line search and model update
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let’s finish off the inversion by running through the line search, which
will generate new models using the gradient, evaluate the objective
function by running forward simulations, and comparing the evaluated
objective function with the value obtained in
**evalaute_initial_misfit**.

Satisfactory reduction in the objective function will result in a
termination of the line search. We are using a bracketing line search
here `(Modrak et
al. 2018) <https://academic.oup.com/gji/article/206/3/1864/2583505>`__,
which requires finding models which both increase and decrease the
misfit with respect to the initial evaluation. Therefore it takes
atleast two trial steps to complete the line search.

.. code:: ipython3

    ! seisflows par stop_after perform_line_search  # We don't want to run the finalize_iteration argument so that we can explore the dir


.. parsed-literal::

    stop_after: evaluate_gradient_from_kernels -> perform_line_search


.. code:: ipython3

    ! seisflows submit


.. parsed-literal::

    2022-08-16 14:41:12 (D) | setting iteration==1 from state file
    2022-08-16 14:41:12 (I) | 
    ================================================================================
                             SETTING UP INVERSION WORKFLOW                          
    ================================================================================
    2022-08-16 14:41:18 (D) | running setup for module 'system.Workstation'
    2022-08-16 14:41:21 (D) | copying par/log file to: /home/bchow/Work/scratch/logs/sflog_003.txt
    2022-08-16 14:41:21 (D) | copying par/log file to: /home/bchow/Work/scratch/logs/parameters_003.yaml
    2022-08-16 14:41:21 (D) | running setup for module 'solver.Specfem2D'
    2022-08-16 14:41:21 (I) | initializing 3 solver directories
    2022-08-16 14:41:22 (D) | running setup for module 'preprocess.Default'
    2022-08-16 14:41:24 (D) | running setup for module 'optimize.Gradient'
    2022-08-16 14:41:26 (I) | re-loading optimization module from checkpoint
    2022-08-16 14:41:28 (I) | re-loading optimization module from checkpoint
    2022-08-16 14:41:28 (I) | 
    ////////////////////////////////////////////////////////////////////////////////
                                  RUNNING ITERATION 01                              
    ////////////////////////////////////////////////////////////////////////////////
    2022-08-16 14:41:28 (I) | 
    ================================================================================
                               RUNNING INVERSION WORKFLOW                           
    ================================================================================
    2022-08-16 14:41:28 (I) | 'evaluate_initial_misfit' has already been run, skipping
    2022-08-16 14:41:28 (I) | 'run_adjoint_simulations' has already been run, skipping
    2022-08-16 14:41:28 (I) | 'postprocess_event_kernels' has already been run, skipping
    2022-08-16 14:41:28 (I) | 'evaluate_gradient_from_kernels' has already been run, skipping
    2022-08-16 14:41:28 (I) | initializing 'bracket'ing line search
    2022-08-16 14:41:28 (I) | enforcing max step length safeguard
    2022-08-16 14:41:28 (D) | step length(s) = 0.00E+00
    2022-08-16 14:41:28 (D) | misfit val(s)  = 1.28E-03
    2022-08-16 14:41:28 (I) | try: first evaluation, attempt guess step length, alpha=9.08E+11
    2022-08-16 14:41:28 (I) | try: applying initial step length safegaurd as alpha has exceeded maximum step length, alpha_new=1.44E+10
    2022-08-16 14:41:28 (D) | overwriting initial step length, alpha_new=2.32E+09
    2022-08-16 14:41:28 (I) | trial model 'm_try' parameters: 
    2022-08-16 14:41:28 (I) | 5800.00 <= vp <= 5800.00
    2022-08-16 14:41:28 (I) | 3244.51 <= vs <= 3790.00
    2022-08-16 14:41:29 (I) | 
    LINE SEARCH STEP COUNT 01
    --------------------------------------------------------------------------------
    2022-08-16 14:41:29 (I) | evaluating objective function for source 001
    2022-08-16 14:41:29 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:41:33 (D) | quantifying misfit with 'Default'
    2022-08-16 14:41:33 (I) | evaluating objective function for source 002
    2022-08-16 14:41:33 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:41:36 (D) | quantifying misfit with 'Default'
    2022-08-16 14:41:36 (I) | evaluating objective function for source 003
    2022-08-16 14:41:36 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:41:40 (D) | quantifying misfit with 'Default'
    2022-08-16 14:41:40 (D) | misfit for trial model (f_try) == 8.65E-04
    2022-08-16 14:41:40 (D) | step length(s) = 0.00E+00, 2.32E+09
    2022-08-16 14:41:40 (D) | misfit val(s)  = 1.28E-03, 8.65E-04
    2022-08-16 14:41:40 (I) | try: misfit not bracketed, increasing step length using golden ratio, alpha=3.76E+09
    2022-08-16 14:41:40 (I) | line search model 'm_try' parameters: 
    2022-08-16 14:41:40 (I) | 5800.00 <= vp <= 5800.00
    2022-08-16 14:41:40 (I) | 3086.61 <= vs <= 3969.23
    2022-08-16 14:41:40 (I) | trial step unsuccessful. re-attempting line search
    2022-08-16 14:41:40 (I) | 
    LINE SEARCH STEP COUNT 02
    --------------------------------------------------------------------------------
    2022-08-16 14:41:40 (I) | evaluating objective function for source 001
    2022-08-16 14:41:40 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:41:44 (D) | quantifying misfit with 'Default'
    2022-08-16 14:41:44 (I) | evaluating objective function for source 002
    2022-08-16 14:41:44 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:41:48 (D) | quantifying misfit with 'Default'
    2022-08-16 14:41:48 (I) | evaluating objective function for source 003
    2022-08-16 14:41:48 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:41:52 (D) | quantifying misfit with 'Default'
    2022-08-16 14:41:52 (D) | misfit for trial model (f_try) == 1.73E-03
    2022-08-16 14:41:52 (D) | step length(s) = 0.00E+00, 2.32E+09, 3.76E+09
    2022-08-16 14:41:52 (D) | misfit val(s)  = 1.28E-03, 8.65E-04, 1.73E-03
    2022-08-16 14:41:52 (I) | try: bracket acceptable but step length unreasonable attempting to re-adjust step length alpha=1.59E+09
    2022-08-16 14:41:52 (I) | line search model 'm_try' parameters: 
    2022-08-16 14:41:52 (I) | 5800.00 <= vp <= 5800.00
    2022-08-16 14:41:52 (I) | 3325.01 <= vs <= 3698.63
    2022-08-16 14:41:52 (I) | trial step unsuccessful. re-attempting line search
    2022-08-16 14:41:52 (I) | 
    LINE SEARCH STEP COUNT 03
    --------------------------------------------------------------------------------
    2022-08-16 14:41:52 (I) | evaluating objective function for source 001
    2022-08-16 14:41:52 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:41:56 (D) | quantifying misfit with 'Default'
    2022-08-16 14:41:56 (I) | evaluating objective function for source 002
    2022-08-16 14:41:56 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:42:00 (D) | quantifying misfit with 'Default'
    2022-08-16 14:42:00 (I) | evaluating objective function for source 003
    2022-08-16 14:42:00 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:42:03 (D) | quantifying misfit with 'Default'
    2022-08-16 14:42:03 (D) | misfit for trial model (f_try) == 2.59E-03
    2022-08-16 14:42:03 (D) | step length(s) = 0.00E+00, 1.59E+09, 2.32E+09, 3.76E+09
    2022-08-16 14:42:03 (D) | misfit val(s)  = 1.28E-03, 2.59E-03, 8.65E-04, 1.73E-03
    2022-08-16 14:42:03 (I) | try: bracket acceptable but step length unreasonable attempting to re-adjust step length alpha=2.82E+09
    2022-08-16 14:42:03 (I) | line search model 'm_try' parameters: 
    2022-08-16 14:42:03 (I) | 5800.00 <= vp <= 5800.00
    2022-08-16 14:42:03 (I) | 3189.77 <= vs <= 3852.13
    2022-08-16 14:42:03 (I) | trial step unsuccessful. re-attempting line search
    2022-08-16 14:42:03 (I) | 
    LINE SEARCH STEP COUNT 04
    --------------------------------------------------------------------------------
    2022-08-16 14:42:03 (I) | evaluating objective function for source 001
    2022-08-16 14:42:03 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:42:07 (D) | quantifying misfit with 'Default'
    2022-08-16 14:42:07 (I) | evaluating objective function for source 002
    2022-08-16 14:42:07 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:42:11 (D) | quantifying misfit with 'Default'
    2022-08-16 14:42:11 (I) | evaluating objective function for source 003
    2022-08-16 14:42:11 (D) | running forward simulation with 'Specfem2D'
    2022-08-16 14:42:15 (D) | quantifying misfit with 'Default'
    2022-08-16 14:42:15 (D) | misfit for trial model (f_try) == 3.46E-03
    2022-08-16 14:42:15 (D) | step length(s) = 0.00E+00, 1.59E+09, 2.32E+09, 2.82E+09, 3.76E+09
    2022-08-16 14:42:15 (D) | misfit val(s)  = 1.28E-03, 2.59E-03, 8.65E-04, 3.46E-03, 1.73E-03
    2022-08-16 14:42:15 (I) | pass: bracket acceptable and step length reasonable. returning minimum line search misfit.
    2022-08-16 14:42:15 (I) | line search model 'm_try' parameters: 
    2022-08-16 14:42:15 (I) | 5800.00 <= vp <= 5800.00
    2022-08-16 14:42:15 (I) | 3244.51 <= vs <= 3790.00
    2022-08-16 14:42:15 (I) | trial step successful. finalizing line search
    2022-08-16 14:42:15 (I) | 
    FINALIZING LINE SEARCH
    --------------------------------------------------------------------------------
    2022-08-16 14:42:15 (I) | writing optimization stats
    2022-08-16 14:42:15 (I) | renaming current (new) optimization vectors as previous model (old)
    2022-08-16 14:42:15 (I) | setting accepted trial model (try) as current model (new)
    2022-08-16 14:42:15 (I) | misfit of accepted trial model is f=8.645E-04
    2022-08-16 14:42:15 (I) | resetting line search step count to 0
    2022-08-16 14:42:15 (I) | stop workflow at `stop_after`: perform_line_search


From the log statements above, we can see that the SeisFlows line search
required 4 trial steps, where it modified values of Vs (shear-wave
velocity) until satisfactory reduction in the objective function was
met. This was the final step in the iteration, and so the finalization
of the line search made preparations for a subsequent iteration.

.. code:: ipython3

    # We can see that we have 'new' and 'old' values for each of the optimization values,
    # representing the previous model (M00) and the current model (M01).
    ! ls scratch/optimize


.. parsed-literal::

    alpha.txt	f_new.txt  f_try.txt  m_new.npz  output_optim.txt
    checkpoint.npz	f_old.txt  g_old.npz  m_old.npz  p_old.npz


.. code:: ipython3

    # The stats/ directory contains text files describing the optimization/line search
    ! cat scratch/optimize/output_optim.txt


.. parsed-literal::

    step_count,step_length,gradient_norm_L1,gradient_norm_L2,misfit,if_restarted,slope,theta
    04,2.323E+09,9.243E-05,1.049E-06,1.279E-03,0,8.263E-13,0.000E+00


4. Conclusions
~~~~~~~~~~~~~~

We’ve now seen how SeisFlows runs an **Inversion** workflow using the
**Specfem2D** solver on a **Workstation** system. More or less, this is
all you need to run SeisFlows with any combination of modules. The
specificities of a system or numerical solver are already handled
internally by SeisFlows, so if you want to use Specmfe3D_Cartesian as
your solver, you would only need to run
``seisflows par solver specfem3d`` at the beginning of your workflow
(you will also need to set up your Specfem3D models, similar to what we
did for Specfem2D here). To run on a slurm system like Chinook
(University of Alaska Fairbanks), you can run
``seisflows par system chinook``.
