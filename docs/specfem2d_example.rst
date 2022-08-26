Specfem2D Workstation Example
=============================

To demonstrate the inversion capabilities of SeisFlows, we will run a
**Specfem2D synthetic-synthetic example** on a **local machine** (tested
on a Linux workstation running CentOS 7, and an Apple Laptop running
macOS 10.14.6). Many of the setup steps here may be unique to our OS and
workstation, but hopefully they may serve as templates for new Users
wanting to explore SeisFlows.

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
    If you do not have a compiled version of SPECFEM2D, then this example will attempt to automatically download and compile SPECFEM2D. This step may fail if you do not have software required by SPECFEM2D, if there are issues with the SPECFEM2D repository itself, or if the configuration and compiling steps fail. If you run any issues, it is recommended that you manually install and compile SPECFEM2D, and directly provide its path to this example problem using the -r or --specfem2d_repo flags (shown below).

Example #1
~~~~~~~~~~

Example #1 runs a 2 iteration inversion using SPECFEM2D, the default
preprocessing module and a gradient descent optimization algorithm.

.. note::
    Example number 1 is meant to **FAIL** during the line search of Iteration #2, after exceeding the maximum allowable line search step count. This is meant to illustrate line search behavior and allow the User to explore a working directory mid-workflow.

.. code:: ipython3

    ! seisflows examples 1  # print example help dialogue


.. parsed-literal::

    1 example: ex1_specfem2d_workstation_inversion
    
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
                                  ///////////////////                               
    This is a [SPECFEM2D] [WORKSTATION] example, which will run an inversion to
    assess misfit between two homogeneous halfspace models with slightly different
    velocities. [3 events, 1 station, 2 iterations]. The inversion is expected to
    fail after the 5th line search step count of the 2nd iteration. The tasks
    involved include:
    
    1. (optional) Download, configure, compile SPECFEM2D
    2. Set up a SPECFEM2D working directory
    3. Generate starting model from Tape2007 example
    4. Generate target model w/ perturbed starting model
    5. Set up a SeisFlows working directory
    6. Run an inversion workflow
    ================================================================================


You can either setup and run the example in separate tasks using the
``examples setup`` and ``submit`` commands. or directly run the example
after setup using the ``examples run`` command. Use the ``-r`` or
``--specfem2d_repo`` flag to point SeisFlows at an existing SPECFEM2D/
repository (with compiled binaries). If not given, SeisFlows will
automatically download, configure and compile SPECFEM2D in your current
working directory.

.. code:: ipython3

    ! seisflows examples setup 1 -r path/to/specfem2d
    ! seisflows submit
    # The above commands are the same as the below
    ! seisflows examples run 1 --specfem2d_repo path/to/specfem2d

A successfully completed example problem will end with the following log
messages:


Example #2
~~~~~~~~~~

Example #2 runs a 1 iteration inversion using SPECFEM2D, the Pyaflowa
preprocessing module and an L-BFGS optimization algorithm. It
successfully completes the line search and is meant to illustrate the
output of the Pyaflowa preprocessing module.

.. code:: ipython3

    ! seisflows examples 2


.. parsed-literal::

    2 example: ex2_specfem2d_workstation_inversion_w_pyatoa
    
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
                                  SEISFLOWS EXAMPLE 2                               
                                  ///////////////////                               
    This is a [SPECFEM2D] [WORKSTATION] example, which will run an inversion to
    assess misfit between a homogeneous halfspace  and checkerboard model using
    Pyatoa for misfit quantification [2 events, 5 stations, 1 iterations]. The tasks
    involved include:
    
    1. (optional) Download, configure, compile SPECFEM2D
    2. Set up a SPECFEM2D working directory
    3. Generate starting model from Tape2007 example
    4. Generate target model w/ perturbed starting model
    5. Set up a SeisFlows working directory
    6. Run an inversion workflow. The line search is expected to attempt 4 evaluations (i01s04)
    ================================================================================


You can run the example with the same command as shown for Example 1:

.. code:: ipython3

    ! seisflows examples run 2 -r path/to/specfem2d
