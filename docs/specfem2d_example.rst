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

.. code:: bash

    ...
    2022-08-25 17:29:16 (I) | 5800.00 <= vp <= 5800.00
    2022-08-25 17:29:16 (I) | 3236.17 <= vs <= 3802.01
    2022-08-25 17:29:16 (I) | trial step unsuccessful. re-attempting line search
    2022-08-25 17:29:16 (I) | 
    LINE SEARCH STEP COUNT 06
    --------------------------------------------------------------------------------
    2022-08-25 17:29:16 (I) | evaluating objective function for source 001
    2022-08-25 17:29:16 (D) | running forward simulation with 'Specfem2D'
    2022-08-25 17:29:20 (D) | quantifying misfit with 'Default'
    2022-08-25 17:29:20 (I) | evaluating objective function for source 002
    2022-08-25 17:29:20 (D) | running forward simulation with 'Specfem2D'
    2022-08-25 17:29:24 (D) | quantifying misfit with 'Default'
    2022-08-25 17:29:24 (I) | evaluating objective function for source 003
    2022-08-25 17:29:24 (D) | running forward simulation with 'Specfem2D'
    2022-08-25 17:29:28 (D) | quantifying misfit with 'Default'
    2022-08-25 17:29:28 (D) | misfit for trial model (f_try) == 7.53E-03
    2022-08-25 17:29:28 (D) | step length(s) = 0.00E+00, 1.47E+08, 2.95E+08, 5.89E+08, 1.18E+09, 2.36E+09, 4.72E+09
    2022-08-25 17:29:28 (D) | misfit val(s)  = 8.65E-04, 7.53E-03, 6.28E-03, 5.02E-03, 3.77E-03, 2.51E-03, 1.26E-03
    2022-08-25 17:29:28 (I) | fail: bracketing line search has failed to reduce the misfit before exceeding `step_count_max`=5
    2022-08-25 17:29:28 (D) | checking gradient/search direction angle, theta:  0.000
    2022-08-25 17:29:28 (C) | 
    ================================================================================
                                   LINE SEARCH FAILED                               
                                   //////////////////                               
    Line search has failed to reduce the misfit and has run out of fallback options.
    Aborting inversion.
    ================================================================================
    EXAMPLE COMPLETED SUCCESFULLY

    
Using the `working directory documentation page <working_directory.html>`__ you can figure out how to navigate around and look at the results of our small inversion problem. We will have a look at a few of the files and directories here. I've run the example problem in a scratch directory but your output directory should look the same.

.. code:: ipython3

    %cd ~/Work/scratch
    ! ls


.. parsed-literal::

    /home/bchow/Work/scratch
    logs	parameters.yaml  sflog.txt    specfem2d
    output	scratch		 sfstate.txt  specfem2d_workdir


In the ``output/`` directory, we can see the updated model from our
first iteration (MODEL_01) and the gradient that was used to create it
(GRADIENT_01). The 2nd iteration produced a gradient (GRADIENT_02), but
was unable to succesfully reduce the misfit during the line search,
which is why we don’t have a MODEL_02.

.. code:: ipython3

    ! ls output
    ! echo
    ! ls output/MODEL_01


.. parsed-literal::

    GRADIENT_01  GRADIENT_02  MODEL_01  MODEL_INIT	MODEL_TRUE
    
    proc000000_vp.bin  proc000000_vs.bin


Because we’re working with SPECFEM2D, we can plot the models and
gradients that were created during our workflow using the
``seisflows plot2d`` command. If we use the ``--savefig`` option we can
also save the output .png files to disk.

.. code:: ipython3

    ! seisflows plot2d GRADIENT_01 vs_kernel --savefig i02_gradient_vs_kernel.png


.. parsed-literal::

    Figure(707.107x707.107)


.. code:: ipython3

    # Because this docs page was made in a Jupyter Notebook, we need to use IPython to open the resulting .png
    from IPython.display import Image
    Image(filename='i02_gradient_vs_kernel.png') 




.. image:: images/specfem2d_example_files/specfem2d_example_15_0.png



Have a look at the `working directory documentation page <working_directory.html>`__ for more detailed explanations of how to navigate the SeisFlows working directory.

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
