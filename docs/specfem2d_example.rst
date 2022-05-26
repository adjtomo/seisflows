Specfem2D workstation example
=============================

To demonstrate the inversion capabilities of SeisFlows, we will run a
**Specfem2D synthetic-synthetic example** on a **local machine** (Linux
workstation running CentOS 7). Many of the setup steps here will likely
be unique to our OS and workstation, but hopefully they may serve as
templates for new Users wanting to explore SeisFlows.

The numerical solver we will use is:
`SPECFEM2D <https://geodynamics.org/cig/software/specfem2d/>`__. We’ll
also be working in our ``seisflows``
`Conda <https://docs.conda.io/en/latest/>`__ environment, see the
installation documentation page for instructions on how to install and
activate the required Conda environment.

--------------

Option 1: Automated run
-----------------------

We have set up this example to run using a single command line argument.
The following command will run an example script which will (1) download
and compile SPECFEM2D, (2) setup a SPECFEM2D working directory to
generate initial and target models, and (3) Run a SeisFlows inversion.

.. warning:: 
    This example attempts to automatically download and compile SPECFEM2D. This step may fail if you are software required by SPECFEM2D, there are issues with the SPECFEM2D repository itself, or the configuration and compiling steps fail. If you run any issues, it is recommended that you manually install and compile SPECFEM2D, and directly provide its path to this example problem when prompted.

.. code:: ipython3

    seisflows examples run 1


.. parsed-literal::

    Run example: ex1_specfem2d_workstation_inversion
    
                                        @@@@@@@@@@                        
                                   .@@@@.    .%&(  %@.          
                                @@@@   @@@@   &@@@@@@ ,%@       
                             @@@@   @@@,  /@@              @    
                            @@@   @@@@   @@@              @     
                          @@@@   @@@@   @@@                @  @ 
                          @@@   @@@@   ,@@@                @ @  
                         @@@@   @@@@    @@@@              @@ @ @
                         @@@@   @@@@@    @@@@@          @@@ @@ @
                         @@@@    @@@@@     @@@@@@@@@@@@@@  @@  @
                          @@@@    @@@@@@        @@@&     @@@  @ 
                          @@@@@     @@@@@@@@         %@@@@#  @@ 
                            @@@@#      @@@@@@@@@@@@@@@@@   @@   
                             &@@@@@          @@@@(       @@&    
                                @@@@@@@             /@@@@       
                                    @@@@@@@@@@@@@@@@@
                                        @@@@@@@@@@          
    
    ================================================================================
                                  SEISFLOWS EXAMPLE 1                              
                                  ////////////////////                              
    This is a [SPECFEM2D] [WORKSTATION] example, which will run 2 iterations of an
    inversion to assess misfit between two homogeneous halfspace models with
    slightly different velocities, 3 sources and 1 receiver. The tasks involved
    include:
    
    1. (optional) Download, configure, compile SPECFEM2D
    2. Set up a SPECFEM2D working directory
    3. Generate starting model from Tape2007 example
    4. Generate target model w/ perturbed starting model
    5. Set up a SeisFlows working directory
    6. Run 2 iterations of an inversion workflow
    ================================================================================
    If you have already downloaded SPECMFE2D, please input its path here. If blank,
    this example will pull the latest version from GitHub and attempt to configure
    and make the binaries: >

--------------

Option 2: Manual run
--------------------

The notebook below details a walkthrough of the automated run shown
above. This is meant for those who want to understand what is going on
under the hood. You are welcome to follow along on your workstation. The
following Table of Contents outlines the steps we will take in this
tutorial:

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

   a. `SF working directory and parameter
      file <#2a.-SF-working-directory-and-parameter-file>`__
   b. `Initialize SF working
      state <#2b.-Initialize-SF-working-state>`__

3. `Run SeisFlows <#2.-Run-SeisFlows>`__

   a. `Forward simulations <#3a.-Forward-simulations>`__
   b. `Exploring the SF directory
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
working directory that will be used by SeisFlows. We’ll use to Python
OS module to do our filesystem processes just to keep everything in
Python, but this can easily be accomplished in bash.

.. code:: ipython3

    import os
    import glob
    import shutil
    import numpy as np

.. code:: ipython3

    # vvv USER MUST EDIT THE FOLLOWING PATHS vvv
    WORKDIR = "/home/bchow/Work/work/sf_specfem2d_example" 
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
long run, and is the SeisFlows dev team’s preferred method of using
SPECFEM.

.. note::
    All SPECFEM2D/3D/3D_GLOBE need to run successfully are the bin/, DATA/, and OUTPUT_FILES/ directories. Everything else in the repository is not mandatory for running binaries.

In this tutorial we will be using the `Tape2007 example
problem <https://github.com/geodynamics/specfem2d/tree/devel/EXAMPLES/Tape2007>`__
to define our **DATA/** directory (last tested 3/9/22, cf893667).

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

    /home/bchow/Work/work/sf_specfem2d_example/specfem2d_workdir
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
     D a t e : 29 - 04 - 2022                                 T i m e  : 12:24:51
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------
    
    see results in directory: OUTPUT_FILES/
    
    done
    Fri Apr 29 12:24:51 AKDT 2022


--------------

Now we need to manually set up our SPECFEM2D working directory. As
mentioned in the previous cell, the only required elements of this
working directory are the following (these files will form the basis for
how SeisFlows operates within the SPECFEM2D framework):

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
    We can use the SeisFlows command line option `seisflows sempar` to directly edit the SPECFEM2D Par_file in the command line. This will work for the SPECFEM3D Par_file as well.

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
    
    # Move the model files (*.bin) into the OUTPUT_FILES directory, where SeisFlows expects them
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
     D a t e : 29 - 04 - 2022                                 T i m e  : 12:25:24
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
     D a t e : 29 - 04 - 2022                                 T i m e  : 12:25:24
     -------------------------------------------------------------------------------
     -------------------------------------------------------------------------------


.. code:: ipython3

    # Great, we have all the necessary SPECFEM files to run our SeisFlows inversion!
    ! ls


.. parsed-literal::

    bin  DATA  OUTPUT_FILES_INIT  OUTPUT_FILES_TRUE


2. Initialize SeisFlows (SF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this Section we will look at a SeisFlows working directory,
parameter file, and working state.

2a. SF working directory and parameter file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As with SPECFEM, SeisFlows requires a parameter file
(**parameters.yaml**) that controls how an automated workflow will
proceed. Because SeisFlows is modular, there are a large number of
potential parameters which may be present in SF parameter file, as each
sub-module may have its own set of unique parameters.

In contrast to SPECFEM’s method of listing all available parameters and
leaving it up the User to determine which ones are relevant to them,
SeisFlows dynamically builds its parameter file based on User inputs.
In this subsection we will use the built-in SeisFlows command line
tools to generate and populate the parameter file.

.. note::
    See the `parameter file documentation page <parameter_file.html>`__ for a more in depth exploration of this central SeisFlows file.

In the previous section we saw the ``sempar`` command in action. We can
use the ``-h`` or help flag to list all available SiesFlows3 command
line commands.

.. code:: ipython3

    ! seisflows -h


.. parsed-literal::

    usage: seisflows [-h] [-w [WORKDIR]] [-p [PARAMETER_FILE]]
                     {setup,configure,init,submit,resume,restart,clean,par,sempar,check,print,convert,reset,debug,edit,examples}
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
        init                Initiate working environment
        submit              Submit initial workflow to system
        resume              Re-submit previous workflow to system
        restart             Remove current environment and submit new workflow
        clean               Remove files relating to an active working environment
        par                 View and edit SeisFlows parameter file
        sempar              View and edit SPECFEM parameter file
        check               Check state of an active environment
        print               Print information related to an active environment
        convert             Convert model file format
        reset               Reset modules within an active state
        debug               Start interactive debug environment
        edit                Open source code file in text editor
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
    parameters.yaml  specfem2d  specfem2d_workdir


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
    #       > seisflows print module
    #   - To auto-fill with docstrings and default values (recommended), run:
    #       > seisflows configure
    #   - To set values as NoneType, use: null
    #   - To set values as infinity, use: inf
    #
    #                                    MODULES
    #                                    ///////
    # WORKFLOW (str):    The method for running SeisFlows; equivalent to main()
    # SOLVER (str):      External numerical solver to use for waveform simulations
    # SYSTEM (str):      Computer architecture of the system being used
    # OPTIMIZE (str):    Optimization algorithm for the inverse problem
    # PREPROCESS (str):  Preprocessing schema for waveform data
    # POSTPROCESS (str): Postprocessing schema for kernels and gradients
    # ==============================================================================
    WORKFLOW: inversion
    SOLVER: specfem2d
    SYSTEM: workstation
    OPTIMIZE: LBFGS 
    PREPROCESS: base
    POSTPROCESS: base


.. code:: ipython3

    # We can use the `seisflows print modules` command to list out the available options 
    ! seisflows print modules


.. parsed-literal::

                                   SEISFLOWS MODULES                               
                                   //////////////////                               
    '+': package, '-': module, '*': class
    
    + SYSTEM
        - seisflows
            * base
            * cluster
            * lsf
            * slurm
            * workstation
        - seisflows-super
            * chinook
            * maui
    + PREPROCESS
        - seisflows
            * base
            * pyatoa
        - seisflows-super
            * pyatoa_nz
    + SOLVER
        - seisflows
            * base
            * specfem2d
            * specfem3d
            * specfem3d_globe
        - seisflows-super
            * specfem3d_maui
    + POSTPROCESS
        - seisflows
            * base
        - seisflows-super
    + OPTIMIZE
        - seisflows
            * LBFGS
            * NLCG
            * base
        - seisflows-super
    + WORKFLOW
        - seisflows
            * base
            * inversion
            * migration
            * test
        - seisflows-super
            * thrifty_inversion
            * thrifty_maui


.. code:: ipython3

    # For this example, we can use most of the default modules, however we need to 
    # change the SOLVER module to let SeisFlows know we're using SPECFEM2D (as opposed to 3D)
    ! seisflows par solver specfem2d
    ! cat parameters.yaml


.. parsed-literal::

    SOLVER: specfem2d -> specfem2d
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
    #       > seisflows print module
    #   - To auto-fill with docstrings and default values (recommended), run:
    #       > seisflows configure
    #   - To set values as NoneType, use: null
    #   - To set values as infinity, use: inf
    #
    #                                    MODULES
    #                                    ///////
    # WORKFLOW (str):    The method for running SeisFlows; equivalent to main()
    # SOLVER (str):      External numerical solver to use for waveform simulations
    # SYSTEM (str):      Computer architecture of the system being used
    # OPTIMIZE (str):    Optimization algorithm for the inverse problem
    # PREPROCESS (str):  Preprocessing schema for waveform data
    # POSTPROCESS (str): Postprocessing schema for kernels and gradients
    # ==============================================================================
    WORKFLOW: inversion
    SOLVER: specfem2d
    SYSTEM: workstation
    OPTIMIZE: LBFGS 
    PREPROCESS: base
    POSTPROCESS: base


--------------

The ``seisflows configure`` command populates the parameter file based
on the chosen modules. SeisFlows will attempt to fill in all parameters
with default values when possible, but values that the User **MUST** set
will be denoted by the value:

   **!!! REQUIRED PARAMETER !!!**

SeisFlows will not work until all of these required parameters are set
by the User. Docstrings above each module show descriptions and
available options for each of these parameters. In the follownig cell we
will use the ``seisflows par`` command to edit the parameters.yaml file
directly, replacing each of the required parameters with a chosen value.
Comments next to each evaluation describe the choice for each.

.. code:: ipython3

    ! seisflows configure
    ! cat parameters.yaml


.. parsed-literal::

    filling parameters.yaml w/ default values
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
    #       > seisflows print module
    #   - To auto-fill with docstrings and default values (recommended), run:
    #       > seisflows configure
    #   - To set values as NoneType, use: null
    #   - To set values as infinity, use: inf
    #
    #                                    MODULES
    #                                    ///////
    # WORKFLOW (str):    The method for running SeisFlows; equivalent to main()
    # SOLVER (str):      External numerical solver to use for waveform simulations
    # SYSTEM (str):      Computer architecture of the system being used
    # OPTIMIZE (str):    Optimization algorithm for the inverse problem
    # PREPROCESS (str):  Preprocessing schema for waveform data
    # POSTPROCESS (str): Postprocessing schema for kernels and gradients
    # ==============================================================================
    WORKFLOW: inversion
    SOLVER: specfem2d
    SYSTEM: workstation
    OPTIMIZE: LBFGS 
    PREPROCESS: base
    POSTPROCESS: base
    
    # =============================================================================
    #                                    SYSTEM                                    
    #                                    //////                                    
    # TITLE (str):
    #   The name used to submit jobs to the system, defaults to the name of the
    #   working directory
    # PRECHECK (list):
    #   A list of parameters that will be displayed to stdout before 'submit' or
    #   'resume' is run. Useful for manually reviewing important parameters prior
    #   to system submission
    # LOG_LEVEL (str):
    #   Verbosity output of SF logger. Available from least to most verbosity:
    #   'CRITICAL', 'WARNING', 'INFO', 'DEBUG'; defaults to 'DEBUG'
    # VERBOSE (bool):
    #   Level of verbosity provided to the output log. If True, log statements
    #   will declare what module/class/function they are being called from.
    #   Useful for debugging but also very noisy.
    # MPIEXEC (str):
    #   Function used to invoke executables on the system. For example 'srun' on
    #   SLURM systems, or './' on a workstation. If left blank, will guess based
    #   on the system.
    # NTASK (int):
    #   Number of separate, individual tasks. Also equal to the number of desired
    #   sources in workflow
    # NPROC (int):
    #   Number of processor to use for each simulation
    # =============================================================================
    TITLE: sf_specfem2d_example
    PRECHECK:
        - TITLE
    LOG_LEVEL: DEBUG
    VERBOSE: False
    MPIEXEC:
    NTASK: 1
    NPROC: 1
    
    # =============================================================================
    #                                  PREPROCESS                                  
    #                                  //////////                                  
    # MISFIT (str):
    #   Misfit function for waveform comparisons, for available see
    #   seisflows.plugins.misfit
    # BACKPROJECT (str):
    #   Backprojection function for migration, for available see
    #   seisflows.plugins.adjoint
    # NORMALIZE (list):
    #   Data normalization parameters used to normalize the amplitudes of
    #   waveforms. Choose from two sets: ENORML1: normalize per event by L1 of
    #   traces; OR ENORML2: normalize per event by L2 of traces; AND TNORML1:
    #   normalize per trace by L1 of itself; OR TNORML2: normalize per trace by
    #   L2 of itself
    # FILTER (str):
    #   Data filtering type, available options are:BANDPASS (req. MIN/MAX
    #   PERIOD/FREQ);LOWPASS (req. MAX_FREQ or MIN_PERIOD); HIGHPASS (req.
    #   MIN_FREQ or MAX_PERIOD)
    # MIN_PERIOD (float):
    #   Minimum filter period applied to time series.See also MIN_FREQ, MAX_FREQ,
    #   if User defines FREQ parameters, they will overwrite PERIOD parameters.
    # MAX_PERIOD (float):
    #   Maximum filter period applied to time series.See also MIN_FREQ, MAX_FREQ,
    #   if User defines FREQ parameters, they will overwrite PERIOD parameters.
    # MIN_FREQ (float):
    #   Maximum filter frequency applied to time series.See also MIN_PERIOD,
    #   MAX_PERIOD, if User defines FREQ parameters, they will overwrite PERIOD
    #   parameters.
    # MAX_FREQ (float):
    #   Maximum filter frequency applied to time series,See also MIN_PERIOD,
    #   MAX_PERIOD, if User defines FREQ parameters, they will overwrite PERIOD
    #   parameters.
    # MUTE (list):
    #   Data mute parameters used to zero out early / late arrivals or offsets.
    #   Choose any number of: EARLY: mute early arrivals; LATE: mute late
    #   arrivals; SHORT: mute short source-receiver distances; LONG: mute long
    #   source-receiver distances
    # =============================================================================
    MISFIT: waveform
    BACKPROJECT: null
    NORMALIZE: []
    FILTER: null
    MIN_PERIOD:
    MAX_PERIOD:
    MIN_FREQ:
    MAX_FREQ:
    MUTE: []
    
    # =============================================================================
    #                                    SOLVER                                    
    #                                    //////                                    
    # MATERIALS (str):
    #   Material parameters used to define model. Available: ['ELASTIC': Vp, Vs,
    #   'ACOUSTIC': Vp, 'ISOTROPIC', 'ANISOTROPIC']
    # DENSITY (str):
    #   How to treat density during inversion. Available: ['CONSTANT': Do not
    #   update density, 'VARIABLE': Update density]
    # ATTENUATION (str):
    #   If True, turn on attenuation during forward simulations, otherwise set
    #   attenuation off. Attenuation is always off for adjoint simulations.
    # COMPONENTS (str):
    #   Components used to generate data, formatted as a single string, e.g. ZNE
    #   or NZ or E
    # SOLVERIO (int):
    #   The format external solver files. Available: ['fortran_binary', 'adios']
    # NT (float):
    #   Number of time steps set in the SPECFEM Par_file
    # DT (float):
    #   Time step or delta set in the SPECFEM Par_file
    # F0 (float):
    #   Dominant source frequency
    # FORMAT (float):
    #   Format of synthetic waveforms used during workflow, available options:
    #   ['ascii', 'su']
    # SOURCE_PREFIX (str):
    #   Prefix of SOURCE files in path SPECFEM_DATA. By default, 'SOURCE' for
    #   SPECFEM2D
    # =============================================================================
    MATERIALS: !!! REQUIRED PARAMETER !!!
    DENSITY: !!! REQUIRED PARAMETER !!!
    ATTENUATION: !!! REQUIRED PARAMETER !!!
    COMPONENTS: ZNE
    SOLVERIO: fortran_binary
    NT: !!! REQUIRED PARAMETER !!!
    DT: !!! REQUIRED PARAMETER !!!
    F0: !!! REQUIRED PARAMETER !!!
    FORMAT: !!! REQUIRED PARAMETER !!!
    SOURCE_PREFIX: SOURCE
    
    # =============================================================================
    #                                  POSTPROCESS                                 
    #                                  ///////////                                 
    # SMOOTH_H (float):
    #   Gaussian half-width for horizontal smoothing in units of meters. If 0.,
    #   no smoothing applied
    # SMOOTH_V (float):
    #   Gaussian half-width for vertical smoothing in units of meters
    # TASKTIME_SMOOTH (int):
    #   Large radii smoothing may take longer than normal tasks. Allocate
    #   additional smoothing task time as a multiple of TASKTIME
    # =============================================================================
    SMOOTH_H: 0.0
    SMOOTH_V: 0.0
    TASKTIME_SMOOTH: 1
    
    # =============================================================================
    #                                   OPTIMIZE                                   
    #                                   ////////                                   
    # LINESEARCH (str):
    #   Algorithm to use for line search, see seisflows.plugins.line_search for
    #   available choices
    # PRECOND (str):
    #   Algorithm to use for preconditioning gradients, see
    #   seisflows.plugins.preconds for available choices
    # STEPCOUNTMAX (int):
    #   Max number of trial steps in line search before a change in line search
    #   behavior
    # STEPLENINIT (float):
    #   Initial line search step length, as a fraction of current model
    #   parameters
    # STEPLENMAX (float):
    #   Max allowable step length, as a fraction of current model parameters
    # LBFGSMEM (int):
    #   Max number of previous gradients to retain in local memory
    # LBFGSMAX (int):
    #   LBFGS periodic restart interval, between 1 and 'inf'
    # LBFGSTHRESH (float):
    #   LBFGS angle restart threshold
    # =============================================================================
    LINESEARCH: Backtrack
    PRECOND:
    STEPCOUNTMAX: 10
    STEPLENINIT: 0.05
    STEPLENMAX: 0.5
    LBFGSMEM: 3
    LBFGSMAX: inf
    LBFGSTHRESH: 0.0
    
    # =============================================================================
    #                                   WORKFLOW                                   
    #                                   ////////                                   
    # CASE (str):
    #   Type of inversion, available: ['data': real data inversion, 'synthetic':
    #   synthetic-synthetic inversion]
    # RESUME_FROM (str):
    #   Name of task to resume inversion from
    # STOP_AFTER (str):
    #   Name of task to stop inversion after finishing
    # SAVEMODEL (bool):
    #   Save final model files after each iteration
    # SAVEGRADIENT (bool):
    #   Save gradient files after each iteration
    # SAVEKERNELS (bool):
    #   Save event kernel files after each iteration
    # SAVETRACES (bool):
    #   Save waveform traces after each iteration
    # SAVERESIDUALS (bool):
    #   Save waveform residuals after each iteration
    # SAVEAS (str):
    #   Format to save models, gradients, kernels. Available: ['binary': save
    #   files in native SPECFEM .bin format, 'vector': save files as NumPy .npy
    #   files, 'both': save as both binary and vectors]
    # BEGIN (int):
    #   First iteration of workflow, 1 <= BEGIN <= inf
    # END (int):
    #   Last iteration of workflow, BEGIN <= END <= inf
    # =============================================================================
    CASE: !!! REQUIRED PARAMETER !!!
    RESUME_FROM:
    STOP_AFTER:
    SAVEMODEL: True
    SAVEGRADIENT: True
    SAVEKERNELS: False
    SAVETRACES: False
    SAVERESIDUALS: False
    SAVEAS: binary
    BEGIN: 1
    END: !!! REQUIRED PARAMETER !!!
    
    # =============================================================================
    #                                     PATHS                                    
    #                                     /////                                    
    # SCRATCH:
    #   scratch path to hold temporary data during workflow
    # OUTPUT:
    #   directory to save workflow outputs to disk
    # SYSTEM:
    #   scratch path to hold any system related data
    # LOCAL:
    #   path to local data to be used during workflow
    # LOGFILE:
    #   the main output log file where all processes will track their status
    # SOLVER:
    #   scratch path to hold solver working directories
    # SPECFEM_BIN:
    #   path to the SPECFEM binary executables
    # SPECFEM_DATA:
    #   path to the SPECFEM DATA/ directory containing the 'Par_file', 'STATIONS'
    #   file and 'CMTSOLUTION' files
    # DATA:
    #   path to data available to workflow
    # MASK:
    #   Directory to mask files for gradient masking
    # OPTIMIZE:
    #   scratch path to store data related to nonlinear optimization
    # MODEL_INIT:
    #   location of the initial model to be used for workflow
    # MODEL_TRUE:
    #   Target model to be used for PAR.CASE == 'synthetic'
    # FUNC:
    #   scratch path to store data related to function evaluations
    # GRAD:
    #   scratch path to store data related to gradient evaluations
    # HESS:
    #   scratch path to store data related to Hessian evaluations
    # =============================================================================
    PATHS:
        SCRATCH: /home/bchow/Work/work/sf_specfem2d_example/scratch
        OUTPUT: /home/bchow/Work/work/sf_specfem2d_example/output
        SYSTEM: /home/bchow/Work/work/sf_specfem2d_example/scratch/system
        LOCAL:
        LOGFILE: /home/bchow/Work/work/sf_specfem2d_example/output_sf.txt
        SOLVER: /home/bchow/Work/work/sf_specfem2d_example/scratch/solver
        SPECFEM_BIN: !!! REQUIRED PATH !!!
        SPECFEM_DATA: !!! REQUIRED PATH !!!
        DATA:
        MASK:
        OPTIMIZE: /home/bchow/Work/work/sf_specfem2d_example/scratch/optimize
        MODEL_INIT: !!! REQUIRED PATH !!!
        MODEL_TRUE:
        FUNC: /home/bchow/Work/work/sf_specfem2d_example/scratch/scratch
        GRAD: /home/bchow/Work/work/sf_specfem2d_example/scratch/evalgrad
        HESS: /home/bchow/Work/work/sf_specfem2d_example/scratch/evalhess


.. code:: ipython3

    # We can check which parameters we will NEED to fill out before running the workflow with the --required flag
    ! seisflows par --required


.. parsed-literal::

    !!! REQUIRED PARAMETER !!!
    ==========================
    	MATERIALS
    	DENSITY
    	ATTENUATION
    	NT
    	DT
    	F0
    	FORMAT
    	CASE
    	END
    !!! REQUIRED PATH !!!
    =====================
    	SPECFEM_BIN
    	SPECFEM_DATA
    	MODEL_INIT


.. code:: ipython3

    # EDIT THE SEISFLOWS PARAMETER FILE
    ! seisflows par materials elastic  # how the velocity model is parameterized
    ! seisflows par density constant  # update density or keep constant
    ! seisflows par attenuation False
    ! seisflows par nt 5000  # set by SPECFEM2D Par_file
    ! seisflows par dt .06  # set by SPECFEM2D Par_file
    ! seisflows par f0 0.084  # set by SOURCE file
    ! seisflows par format ascii  # how to output synthetic seismograms
    ! seisflows par begin 1  # first iteration
    ! seisflows par end 1  # final iteration -- we will only run 1
    ! seisflows par case synthetic  # synthetic-synthetic means we need both INIT and TRUE models
    
    # Use Python syntax here to access path constants
    os.system(f"seisflows par specfem_bin {SPECFEM2D_BIN}")  # set path to SPECFEM2D binaries
    os.system(f"seisflows par specfem_data {SPECFEM2D_DATA}")  # set path to SEPCFEM2D DATA/
    os.system(f"seisflows par model_init {SPECFEM2D_MODEL_INIT}")  # set path to INIT model
    os.system(f"seisflows par model_true {SPECFEM2D_MODEL_TRUE}")  # set path to TRUE model


.. parsed-literal::

    MATERIALS: !!! REQUIRED PARAMETER !!! -> elastic
    DENSITY: !!! REQUIRED PARAMETER !!! -> constant
    ATTENUATION: !!! REQUIRED PARAMETER !!! -> False
    NT: !!! REQUIRED PARAMETER !!! -> 5000
    DT: !!! REQUIRED PARAMETER !!! -> .06
    F0: !!! REQUIRED PARAMETER !!! -> 0.084
    FORMAT: !!! REQUIRED PARAMETER !!! -> ascii
    BEGIN: 1 -> 1
    END: !!! REQUIRED PARAMETER !!! -> 1
    CASE: !!! REQUIRED PARAMETER !!! -> synthetic




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


2b. Initialize SF working state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The SeisFlows command ``seisflows init`` will generate the a SeisFlows
working state without submitting any jobs to the system. This is useful
for testing to see if the user has set an acceptable parameter file, and
if SeisFlows is working as expected.

The result of running ``seisflows init`` is a collection of pickle (*.p)
and JSON files which define the active Python environment. SeisFlows
relies directly on these files to determine where it is in a workflow.
Throughout an active workflow, SeisFlows will checkpoint itself to
these pickle and JSON files such that if a workflow finishes or crashes,
the User can resume a workflow from the last checkpointed state rather
than needing to restart the workflow.

   **DEBUG MODE:** After running ``seisflows init`` you can explore the
   SeisFlows working state in an interactive iPython environment by
   running ``seisflows debug``. This will open up an iPython environment
   in which the active working state is loaded and accessible The debug
   mode is invaluable for exploring the SeisFlows working state,
   debugging errors, and performing manual manipulations to an otherwise
   automated tool. You can try for yourself by running debug mode and
   typing ‘preprocess’ to access the active preprocess module.

.. code:: ipython3

    os.chdir(WORKDIR)
    ! seisflows init
    ! ls output


.. parsed-literal::

    instantiating SeisFlows working state in directory: output
    seisflows_optimize.p	   seisflows_postprocess.p  seisflows_system.p
    seisflows_parameters.json  seisflows_preprocess.p   seisflows_workflow.p
    seisflows_paths.json	   seisflows_solver.p


.. code:: ipython3

    # All of the parameters defined in parameters.yaml are saved in this 
    # internally-used JSON file
    ! head output/seisflows_parameters.json


.. parsed-literal::

    {
        "ATTENUATION": false,
        "BACKPROJECT": null,
        "BEGIN": 1,
        "CASE": "synthetic",
        "COMPONENTS": "ZNE",
        "DENSITY": "constant",
        "DT": 0.06,
        "END": 1,
        "F0": 0.084,


.. code:: ipython3

    # Similarly, paths that SeisFlows uses to navigate the system are stored
    # in the seisflows_paths.json file
    ! head output/seisflows_paths.json


.. parsed-literal::

    {
        "DATA": null,
        "FUNC": "/home/bchow/Work/work/sf_specfem2d_example/scratch/scratch",
        "GRAD": "/home/bchow/Work/work/sf_specfem2d_example/scratch/evalgrad",
        "HESS": "/home/bchow/Work/work/sf_specfem2d_example/scratch/evalhess",
        "LOCAL": null,
        "LOGFILE": "/home/bchow/Work/work/sf_specfem2d_example/output_sf.txt",
        "MASK": null,
        "MODEL_INIT": "/home/bchow/Work/work/sf_specfem2d_example/specfem2d_workdir/OUTPUT_FILES_INIT",
        "MODEL_TRUE": "/home/bchow/Work/work/sf_specfem2d_example/specfem2d_workdir/OUTPUT_FILES_TRUE",


3. Run SeisFlows
~~~~~~~~~~~~~~~~~

In this Section we will run SeisFlows to generate synthetic
seismograms, kernels, a gradient, and an updated velocity model.

3a. Forward simulations
^^^^^^^^^^^^^^^^^^^^^^^

SeisFlows is an automated workflow tool, such that once we run
``seisflows submit`` we should not need to intervene in the workflow.
However the package does allow the User flexibility in how they want the
workflow to behave.

For example, we can run our workflow in stages by taking advantage of
the ``stop_after`` and ``resume_from`` parameters. As their names
suggest, these parameters allow us to stop and resume the workflow at
certain stages (i.e., functions in workflow.main()).

The available arguments for ``stop_after`` and ``resume_from`` are
discovered by running the command: ``seisflows print flow``, which tells
us what functions will be run from main().

.. code:: ipython3

    ! seisflows print flow


.. parsed-literal::

                                SEISFLOWS WORKFLOW MAIN                            
                                ////////////////////////                            
    Flow arguments for <class 'seisflows.workflow.inversion.Inversion'>
    
    1: setup
    2: initialize
    3: evaluate_gradient
    4: write_gradient
    5: compute_direction
    6: line_search
    7: finalize
    8: clean


--------------

In an inversion (the workflow we have selected) the flow arguments are
described as:

0. **setup:** Not technically listed in the flow arguments, runs setup()
   for all SeisFlows modules. If running a synthetic-synthetic
   workflow, solver.setup() will generate “data” by running the forward
   solver using MODEL_TRUE
1. **initialize:**

   a. Call numerical solver to run forward simulations using MODEL_INIT,
      generating synthetics
   b. Evaluate the objective function by performing waveform comparisons
   c. Prepare ``evaluate gradient`` step by generating adjoint sources
      and auxiliary files

2. **evaluate_gradient:** Call numerical solver to run adjoint
   simulation, generating kernels
3. **write_gradient:** Combine all event kernels into a misfit kernel.
   Optionally smooth and mask the misfit kernel
4. **compute_direction:** Call on the optimization library to scale the
   misfit kernel into the gradient and compute a search direction
5. **line_search:** Perform a line search by algorithmically scaling the
   gradient and evaluating the misfit function (forward simulations and
   misfit quantification) until misfit is acceptably reduced
6. **finalize:** Run any finalization steps such as saving traces,
   kernels, gradients and models to disk, setting up SeisFlows for any
   subsequent iterations.
7. **clean:** Clean the scratch/ directory in preparation for subsequent
   i

Let’s set the ``stop_after`` argument to **initialize**, this will halt
the workflow after the intialization step. We’ll also set the
``verbose`` parameter to ‘False’, to keep the logging format relatively
simple. We will explore the ``verbose``\ ==True option in a later cell.

.. code:: ipython3

    ! seisflows par stop_after initialize
    ! seisflows par verbose False


.. parsed-literal::

    STOP_AFTER:  -> initialize
    VERBOSE: False -> False


--------------

Now let’s run SeisFlows. There are a few ways to do this: ``submit``,
``resume``, and ``restart``

1. Since we already ran ``seisflows init``, the ``seisflows submit``
   option will not work, as SeisFlows considers this an active working
   state and ``submit`` can only be run on uninitialized working states.
2. To run a workflow in an active working state ``resume`` will load the
   current working state from the output/ directory and submit a
   workflow given the current parameter file.
3. The ``restart`` command is simply a convenience function that runs
   ``clean`` (to remove an active working state) and ``submit`` (to
   submit a fresh working state).

Since we haven’t done anything in this working state, we will go with a
modified version of Option 3 by running ``clean`` and then ``submit``.
We’ll use the ``-f`` flag (stands for **‘force’**) to skip over the
standard input prompt that asks the User if they are sure they want to
clean and submit.

But first we’ll try to run ``seisflows submit`` to show why Option 1
**will not work**.

.. code:: ipython3

    ! seisflows submit -f


.. parsed-literal::

    2022-04-29 12:32:17 | initializing SeisFlows in sys.modules
    ================================================================================
                                        WARNING                                     
                                        ///////                                     
    Data from previous workflow found in working directory.
    
    > seisflows restart: delete data and start new workflow
    > seisflows resume: resume existing workflow
    ================================================================================


--------------

**Okay, let’s go!** In the following cell we will run the SeisFlows
Inversion workflow. In the output cell we will see the logging
statements outputted by SeisFlows, both to stdout and to the output log
file (defaults to ./output_seisflows.txt) which details the progress of
our inversion

.. code:: ipython3

    ! seisflows clean -f
    ! seisflows submit -f


.. parsed-literal::

    ================================================================================
                                         CLEAN                                      
                                         /////                                      
    + skipping over: /home/bchow/Work/work/sf_specfem2d_example/parameters.yaml
    - deleting file/folder: /home/bchow/Work/work/sf_specfem2d_example/scratch
    - deleting file/folder: /home/bchow/Work/work/sf_specfem2d_example/stats
    - deleting file/folder: /home/bchow/Work/work/sf_specfem2d_example/output
    - deleting file/folder: /home/bchow/Work/work/sf_specfem2d_example/output_sf.txt
    - deleting file/folder: /home/bchow/Work/work/sf_specfem2d_example/logs
    ================================================================================
    2022-04-29 12:38:37 | initializing SeisFlows in sys.modules
    2022-04-29 12:38:42 | copying par/log file to: /home/bchow/Work/work/sf_specfem2d_example/logs/output_sf_001.txt
    2022-04-29 12:38:42 | copying par/log file to: /home/bchow/Work/work/sf_specfem2d_example/logs/parameters_001.yaml
    2022-04-29 12:38:42 | exporting current working environment to disk
    2022-04-29 12:38:42 | 
    ////////////////////////////////////////////////////////////////////////////////
                      WORKFLOW WILL STOP AFTER FUNC: 'initialize'                   
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:38:42 | 
    ================================================================================
                              STARTING INVERSION WORKFLOW                           
    ================================================================================
    2022-04-29 12:38:42 | 
    ////////////////////////////////////////////////////////////////////////////////
                                    ITERATION 1 / 1                                 
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:38:42 | 
    ////////////////////////////////////////////////////////////////////////////////
                                PERFORMING MODULE SETUP                             
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:38:42 | misfit function is: 'waveform'
    2022-04-29 12:38:43 | writing line search history file:
    /home/bchow/Work/work/sf_specfem2d_example/stats/line_search.txt
    2022-04-29 12:38:44 | checking poissons ratio for: 'm_new.npy'
    2022-04-29 12:38:44 | model parameters (m_new.npy i01s00):
    2022-04-29 12:38:44 | 5800.00 <= vp <= 5800.00
    2022-04-29 12:38:44 | 3500.00 <= vs <= 3500.00
    2022-04-29 12:38:44 | 0.21 <= pr <= 0.21
    2022-04-29 12:38:46 | setting up solver on system...
    2022-04-29 12:38:46 | checkpointing working environment to disk
    2022-04-29 12:38:47 | exporting current working environment to disk
    2022-04-29 12:38:48 | running task solver.setup 1 times
    2022-04-29 12:38:48 | initializing 1 solver directories
    2022-04-29 12:38:53 | source 001 symlinked as mainsolver
    2022-04-29 12:38:53 | generating 'data' with MODEL_TRUE synthetics
    2022-04-29 12:39:00 | running mesh generation for MODEL_INIT
    2022-04-29 12:39:02 | 
    ================================================================================
                                 INITIALIZING INVERSION                             
    ================================================================================
    2022-04-29 12:39:02 | 
    EVALUATE OBJECTIVE FUNCTION
    --------------------------------------------------------------------------------
    2022-04-29 12:39:02 | saving model 'm_new.npy' to:
    /home/bchow/Work/work/sf_specfem2d_example/scratch/evalgrad/model
    2022-04-29 12:39:03 | evaluating objective function 1 times on system...
    2022-04-29 12:39:03 | checkpointing working environment to disk
    2022-04-29 12:39:05 | exporting current working environment to disk
    2022-04-29 12:39:05 | running task solver.eval_func 1 times
    2022-04-29 12:39:05 | running forward simulations
    2022-04-29 12:39:11 | calling preprocess.prepare_eval_grad()
    2022-04-29 12:39:11 | preparing files for gradient evaluation
    2022-04-29 12:39:11 | exporting residuals to:
    /home/bchow/Work/work/sf_specfem2d_example/scratch/evalgrad
    2022-04-29 12:39:12 | summing residuals with preprocess module
    2022-04-29 12:39:12 | saving misfit 1.748E-03 to tag 'f_new.txt'
    2022-04-29 12:39:12 | 
    ================================================================================
                                FINISHED FLOW EXECUTION                             
    ================================================================================
    2022-04-29 12:39:12 | 
    ================================================================================
                              FINISHED INVERSION WORKFLOW                           
    ================================================================================


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


.. code:: ipython3

    # We can also see that we have generated a STATIONS_ADJOINT file, which is required for 
    # running the adjoint simulations (i.e., evaluate the gradient)
    ! head scratch/solver/001/DATA/STATIONS_ADJOINT


.. parsed-literal::

    S0001    AA       180081.4100000       388768.7100000       0.0         0.0


3b. Adjoint simulations
^^^^^^^^^^^^^^^^^^^^^^^

Now that we have all the required files for running an adjoint
simulation (*.adj waveforms and STATIONS_ADJOINT file), we can continue
with the SeisFlows Inversion workflow. No need to edit the Par_file or
anything like that, SeisFlows will take care of that under the hood. We
simply need to tell the workflow (via the parameters.yaml file) to
``resume_from`` the correct function. We can have a look at these
functions again:

.. code:: ipython3

    ! seisflows print flow


.. parsed-literal::

                                SEISFLOWS WORKFLOW MAIN                            
                                ////////////////////////                            
    Flow arguments for <class 'seisflows.workflow.inversion.Inversion'>
    
    1: setup
    2: initialize
    3: evaluate_gradient
    4: write_gradient
    5: compute_direction
    6: line_search
    7: finalize
    8: clean


.. code:: ipython3

    # We'll stop just before the line search so that we can take a look at the files 
    # generated during the middle tasks
    ! seisflows par resume_from evaluate_gradient
    ! seisflows par stop_after compute_direction


.. parsed-literal::

    RESUME_FROM:  -> evaluate_gradient
    STOP_AFTER: initialize -> compute_direction


.. code:: ipython3

    # We can use the `seisflows resume` command to continue an active workflow
    # again we use the '-f' flag to skip past the user-input stage.
    ! seisflows resume -f


.. parsed-literal::

    2022-04-29 12:41:21 | copying par/log file to: /home/bchow/Work/work/sf_specfem2d_example/logs/output_sf_002.txt
    2022-04-29 12:41:21 | copying par/log file to: /home/bchow/Work/work/sf_specfem2d_example/logs/parameters_002.yaml
    2022-04-29 12:41:21 | exporting current working environment to disk
    2022-04-29 12:41:21 | 
    ////////////////////////////////////////////////////////////////////////////////
                  WORKFLOW WILL RESUME FROM FUNC: 'evaluate_gradient'               
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:41:21 | 
    ////////////////////////////////////////////////////////////////////////////////
                   WORKFLOW WILL STOP AFTER FUNC: 'compute_direction'               
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:41:21 | 
    ================================================================================
                              STARTING INVERSION WORKFLOW                           
    ================================================================================
    2022-04-29 12:41:21 | 
    ////////////////////////////////////////////////////////////////////////////////
                                    ITERATION 1 / 1                                 
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:41:21 | 
    ////////////////////////////////////////////////////////////////////////////////
                                  EVALUATING GRADIENT                               
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:41:21 | evaluating gradient 1 times on system...
    2022-04-29 12:41:21 | checkpointing working environment to disk
    2022-04-29 12:41:22 | exporting current working environment to disk
    2022-04-29 12:41:23 | running task solver.eval_grad 1 times
    2022-04-29 12:41:23 | running adjoint simulations
    2022-04-29 12:41:38 | exporting kernels to:
    /home/bchow/Work/work/sf_specfem2d_example/scratch/evalgrad
    2022-04-29 12:41:38 | 
    ////////////////////////////////////////////////////////////////////////////////
                                 POSTPROCESSING KERNELS                             
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:41:38 | processing kernels into gradient on system...
    2022-04-29 12:41:38 | checkpointing working environment to disk
    2022-04-29 12:41:39 | exporting current working environment to disk
    2022-04-29 12:41:39 | running task postprocess.process_kernels 1 times
    2022-04-29 12:41:39 | saving summed kernels to:
    /home/bchow/Work/work/sf_specfem2d_example/scratch/evalgrad/kernels/sum
    2022-04-29 12:41:41 | 
    ////////////////////////////////////////////////////////////////////////////////
                               COMPUTING SEARCH DIRECTION                           
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:41:41 | computing search direction with L-BFGS
    2022-04-29 12:41:41 | first L-BFGS iteration, setting search direction as inverse gradient
    2022-04-29 12:41:41 | 
    ================================================================================
                                FINISHED FLOW EXECUTION                             
    ================================================================================
    2022-04-29 12:41:41 | 
    ================================================================================
                              FINISHED INVERSION WORKFLOW                           
    ================================================================================


--------------

The functions **evaluate_gradient()** through **compute_direction()**
have run adjoint simulations to generate event kernels and sum the
kernels into the misfit kernel.

   **NOTE**: Because we only have one event, our misfit kernel is just
   exactly our event kernel. And since we did not specify any smoothing
   lenghts (PAR.SMOOTH_H and PAR.SMOOTH_V), no smoothing of the gradient
   has occurred.

Using the L-BFGS optimization algorithm, SeisFlows has computed a
search direction that will be used in the line search to search for a
best fitting model which optimally reduces the objective function. We
can take a look at where SeisFlows has stored the information relating
to kernel generation and the optimization computation.

.. code:: ipython3

    # Gradient evaluation files are stored here, the kernels are stored separately from the gradient incase
    # the user wants to manually manipulate them
    ! ls scratch/evalgrad


.. parsed-literal::

    gradient  kernels  model  residuals


.. code:: ipython3

    # SeisFlows stores all kernels and gradient information as SPECFEM binary (.bin) files
    ! ls scratch/evalgrad/gradient


.. parsed-literal::

    proc000000_vp_kernel.bin  proc000000_vs_kernel.bin


.. code:: ipython3

    # Kernels are stored on a per-event basis, and summed together (sum/). If smoothing was performed, 
    # we would see both smoothed and unsmoothed versions of the misfit kernel
    ! ls scratch/evalgrad/kernels


.. parsed-literal::

    001  sum


.. code:: ipython3

    # We can see that some new values have been stored in prepartion for the line search,
    # including g_new (current gradient) and p_new (current search direction). These are also
    # stored as vector NumPy arrays (.npy files)
    ! ls scratch/optimize


.. parsed-literal::

    f_new.txt  g_new.npy  LBFGS  m_new.npy	p_new.npy


.. code:: ipython3

    p_new = np.load("scratch/optimize/p_new.npy")
    print(p_new)


.. parsed-literal::

    [-0.00000000e+00 -0.00000000e+00 -0.00000000e+00 ... -3.96447909e-11
     -2.00156454e-11 -2.61676726e-12]


--------------

3c. Line search and model update
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let’s finish off the inversion by running through the line search, which
will generate new models using the gradient, evaluate the objective
function by running forward simulations, and comparing the evaluated
objective function with the value obtained in **initialize**.
Satisfactory reduction in the objective function will result in a
termination of the line search. We are using a bracketing line search
here (CITE RYANS PAPER), which requires finding models which both
increase and decrease the misfit with respect to the initial evaluation.
Therefore it will likely take more than two trial steps to complete the
line search

.. code:: ipython3

    ! seisflows par resume_from line_search  # resume from the line search 
    ! seisflows par stop_after finalize  # We don't want to run the clean() argument so that we can explore the dir


.. parsed-literal::

    RESUME_FROM: evaluate_gradient -> line_search
    STOP_AFTER: compute_direction -> finalize


.. code:: ipython3

    ! seisflows resume -f


.. parsed-literal::

    2022-04-29 12:42:40 | copying par/log file to: /home/bchow/Work/work/sf_specfem2d_example/logs/output_sf_003.txt
    2022-04-29 12:42:40 | copying par/log file to: /home/bchow/Work/work/sf_specfem2d_example/logs/parameters_003.yaml
    2022-04-29 12:42:40 | exporting current working environment to disk
    2022-04-29 12:42:42 | 
    ////////////////////////////////////////////////////////////////////////////////
                     WORKFLOW WILL RESUME FROM FUNC: 'line_search'                  
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:42:42 | 
    ////////////////////////////////////////////////////////////////////////////////
                       WORKFLOW WILL STOP AFTER FUNC: 'finalize'                    
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:42:42 | 
    ================================================================================
                              STARTING INVERSION WORKFLOW                           
    ================================================================================
    2022-04-29 12:42:42 | 
    ////////////////////////////////////////////////////////////////////////////////
                                    ITERATION 1 / 1                                 
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:42:42 | 
    ================================================================================
                            CONDUCTING LINE SEARCH (i01s00)                         
    ================================================================================
    2022-04-29 12:42:42 | max step length safeguard is: 5.26E+10
    2022-04-29 12:42:42 | 
    EVALUATE BRACKETING LINE SEARCH
    --------------------------------------------------------------------------------
    2022-04-29 12:42:42 | step length(s) = 0.00E+00
    2022-04-29 12:42:42 | misfit val(s)  = 1.75E-03
    2022-04-29 12:42:42 | first iteration, guessing trial step
    2022-04-29 12:42:42 | initial step length safegaurd, setting manual step length
    2022-04-29 12:42:42 | manually set initial step length: 5.26E+09
    2022-04-29 12:42:42 | checking poissons ratio for: 'm_try.npy'
    2022-04-29 12:42:42 | model parameters (m_try.npy i01s00):
    2022-04-29 12:42:42 | 5800.00 <= vp <= 5800.00
    2022-04-29 12:42:42 | 3278.69 <= vs <= 3790.00
    2022-04-29 12:42:42 | 0.13 <= pr <= 0.27
    2022-04-29 12:42:42 | 
    ////////////////////////////////////////////////////////////////////////////////
                                TRIAL STEP COUNT: i01s01                            
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:42:42 | 
    EVALUATE OBJECTIVE FUNCTION
    --------------------------------------------------------------------------------
    2022-04-29 12:42:42 | saving model 'm_try.npy' to:
    /home/bchow/Work/work/sf_specfem2d_example/scratch/scratch/model
    2022-04-29 12:42:42 | evaluating objective function 1 times on system...
    2022-04-29 12:42:42 | checkpointing working environment to disk
    2022-04-29 12:42:44 | exporting current working environment to disk
    2022-04-29 12:42:44 | running task solver.eval_func 1 times
    2022-04-29 12:42:44 | running forward simulations
    2022-04-29 12:42:49 | calling preprocess.prepare_eval_grad()
    2022-04-29 12:42:49 | preparing files for gradient evaluation
    2022-04-29 12:42:50 | exporting residuals to:
    /home/bchow/Work/work/sf_specfem2d_example/scratch/scratch
    2022-04-29 12:42:50 | summing residuals with preprocess module
    2022-04-29 12:42:50 | saving misfit 9.850E-04 to tag 'f_try.txt'
    2022-04-29 12:42:50 | 
    EVALUATE BRACKETING LINE SEARCH
    --------------------------------------------------------------------------------
    2022-04-29 12:42:50 | step length(s) = 0.00E+00, 5.26E+09
    2022-04-29 12:42:50 | misfit val(s)  = 1.75E-03, 9.85E-04
    2022-04-29 12:42:50 | misfit not bracketed, increasing step length
    2022-04-29 12:42:50 | checking poissons ratio for: 'm_try.npy'
    2022-04-29 12:42:50 | model parameters (m_try.npy i01s01):
    2022-04-29 12:42:50 | 5800.00 <= vp <= 5800.00
    2022-04-29 12:42:50 | 3141.92 <= vs <= 3969.23
    2022-04-29 12:42:50 | 0.06 <= pr <= 0.29
    2022-04-29 12:42:50 | retrying with new trial step
    2022-04-29 12:42:50 | 
    ////////////////////////////////////////////////////////////////////////////////
                                TRIAL STEP COUNT: i01s02                            
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 12:42:50 | 
    EVALUATE OBJECTIVE FUNCTION
    --------------------------------------------------------------------------------
    2022-04-29 12:42:50 | saving model 'm_try.npy' to:
    /home/bchow/Work/work/sf_specfem2d_example/scratch/scratch/model
    2022-04-29 12:42:51 | evaluating objective function 1 times on system...
    2022-04-29 12:42:51 | checkpointing working environment to disk
    2022-04-29 12:42:52 | exporting current working environment to disk
    2022-04-29 12:42:53 | running task solver.eval_func 1 times
    2022-04-29 12:42:53 | running forward simulations
    2022-04-29 12:42:59 | calling preprocess.prepare_eval_grad()
    2022-04-29 12:42:59 | preparing files for gradient evaluation
    2022-04-29 12:42:59 | exporting residuals to:
    /home/bchow/Work/work/sf_specfem2d_example/scratch/scratch
    2022-04-29 12:43:00 | summing residuals with preprocess module
    2022-04-29 12:43:00 | saving misfit 1.227E-03 to tag 'f_try.txt'
    2022-04-29 12:43:00 | 
    EVALUATE BRACKETING LINE SEARCH
    --------------------------------------------------------------------------------
    2022-04-29 12:43:00 | step length(s) = 0.00E+00, 5.26E+09, 8.51E+09
    2022-04-29 12:43:00 | misfit val(s)  = 1.75E-03, 9.85E-04, 1.23E-03
    2022-04-29 12:43:00 | bracket okay, step length reasonable, pass
    2022-04-29 12:43:00 | checking poissons ratio for: 'm_try.npy'
    2022-04-29 12:43:00 | model parameters (m_try.npy i01s02):
    2022-04-29 12:43:00 | 5800.00 <= vp <= 5800.00
    2022-04-29 12:43:00 | 3278.69 <= vs <= 3790.00
    2022-04-29 12:43:00 | 0.13 <= pr <= 0.27
    2022-04-29 12:43:00 | trial step successful
    2022-04-29 12:43:00 | 
    FINALIZING LINE SEARCH
    --------------------------------------------------------------------------------
    2022-04-29 12:43:00 | shifting current model (new) to previous model (old)
    2022-04-29 12:43:00 | setting accepted line search model as current model
    2022-04-29 12:43:00 | current misfit is f_new.txt=9.850E-04
    2022-04-29 12:43:00 | writing optimization stats to: stats
    2022-04-29 12:43:00 | resetting line search step count to 0
    2022-04-29 12:43:00 | 
    ================================================================================
                                 FINALIZING ITERATION 1                             
    ================================================================================
    2022-04-29 12:43:00 | exporting current working environment to disk
    2022-04-29 12:43:01 | saving model 'm_new.npy' to path:
    /home/bchow/Work/work/sf_specfem2d_example/output/model_0001
    2022-04-29 12:43:02 | saving gradient to path:
    /home/bchow/Work/work/sf_specfem2d_example/output/gradient_0001
    2022-04-29 12:43:02 | 
    ================================================================================
                                FINISHED FLOW EXECUTION                             
    ================================================================================
    2022-04-29 12:43:02 | 
    ================================================================================
                              FINISHED INVERSION WORKFLOW                           
    ================================================================================


From the log statements above, we can see that the SeisFlows line
search required 2 trial steps, where it modified values of Vs until
satisfactory reduction in the objective function was met. This was the
final step in the iteration, and so the finalization step made
last-minute preparations for a subsequent iteration.

.. code:: ipython3

    # We can see that we have 'new' and 'old' values for each of the optimization values,
    # representing the previous model (M00) and the current model (M01).
    ! ls scratch/optimize


.. parsed-literal::

    alpha.npy  f_old.txt  g_old.npy  m_new.npy  p_old.npy
    f_new.txt  f_try.txt  LBFGS	 m_old.npy


.. code:: ipython3

    # The stats/ directory contains text files describing the optimization/line search
    ! ls stats


.. parsed-literal::

    factor.txt	      line_search.txt  slope.txt	theta.txt
    gradient_norm_L1.txt  misfit.txt       step_count.txt
    gradient_norm_L2.txt  restarted.txt    step_length.txt


.. code:: ipython3

    # For example we can look at the step length chosen for the accepted trial step in the line search
    ! cat stats/line_search.txt


.. parsed-literal::

          ITER     STEPLEN      MISFIT
    ==========  ==========  ==========
             1   0.000e+00   1.748e-03
                 5.261e+09   9.850e-04
                 8.512e+09   1.227e-03


4. Conclusions
~~~~~~~~~~~~~~

We’ve now seen how SeisFlows runs an **Inversion** workflow using the
**Specfem2D** solver on a **serial** system (local workstation). More or
less, this is all you need to run SeisFlows with any combination of
modules. The specificities of a system or numerical solver are already
handled internally by SeisFlows, so if you want to use
Specmfe3D_Cartesian as your solver, you would only need to run
``seisflows par solver specfem3d`` at the beginning of your workflow
(you will also need to setup your Specfem3D models, similar to what we
did for Specfem2D here). To run on a slurm system like Chinook, you can
run ``seisflows par system chinook``.

