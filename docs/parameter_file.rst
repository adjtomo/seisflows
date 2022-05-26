Parameter File
==============

The parameter file is the central control object for a SeisFlows
workflow. Here we take a look at the anatomy of a parameter file.
Parameter files in SeisFlows are formatted in the `YAML format (YAML
Ain’t Markup Language) <https://pyyaml.org/wiki/PyYAMLDocumentation>`__.

Template
--------

Each workflow starts with the module-only template parameter file which
defines the core modules which make up the package. Running
``seisflows setup`` from the command line will create this file.

.. code:: ipython3

    ! seisflows setup


.. parsed-literal::

    creating parameter file: parameters.yaml


.. code:: ipython3

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


How do I choose my modules?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

As seen above, each of the modules comes with a default value. But you
may want to run a migration, not an inversion. Or run with SPECFEM3D not
2D. As stated in the comments at the top of the file, the
``seisflows print modules`` command lists out all available options.
Don’t see an option that works for you? Learn to extend the SeisFlows
package here: **!!! docs page link here !!!**

.. code:: ipython3

    ! seisflows print modules


.. parsed-literal::

                                   SEISFLOWS3 MODULES                               
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


How do I change modules?
~~~~~~~~~~~~~~~~~~~~~~~~

Feel free to use any old text editor to edit the YAML file, or you can
use the ``seisflows par`` command to make changes directly from the
command line. For example, say we want to use SPECFEM3D

.. code:: ipython3

    # Changes the current parameter to the given value
    ! seisflows par solver specfem3d


.. parsed-literal::

    SOLVER: specfem2d -> specfem3d


.. code:: ipython3

    # Prints out the current parameter value
    ! seisflows par solver


.. parsed-literal::

    SOLVER: specfem3d


How do I get to a full parameter file?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module-only parameter file serves as as a template for dynamically
generating a full parameter file. Since each module requires it’s own
unique set of parameters and paths, each parameter file will look
different. We can use the ``seisflows configure`` command to complete
our parmater file, based on the chosen modules.

.. code:: ipython3

    ! seisflows configure


.. parsed-literal::

    filling parameters.yaml w/ default values


Anatomy of the parameter file
-----------------------------

As we will see below, the parameter file has now been generated. Each
module will define its own section, separated by a header of comments.
Within each header, parameter names, types and descriptions are listed.
At the bottom of the parameter file, there is a section defining paths
required by the workflow. Section headers will look something:

.. code:: ipython3

    # =============================================================================
    #                                    MODULE
    #                                    //////                                    
    # PARAMETER_NAME (type):
    #   Description
    # ...
    # =============================================================================
    PARAMETER_NAME: parameter_value

.. code:: ipython3

    ! head -80 parameters.yaml


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
    # WORKFLOW (str):    The method for running SeisFlows; equivalent to main()
    # SOLVER (str):      External numerical solver to use for waveform simulations
    # SYSTEM (str):      Computer architecture of the system being used
    # OPTIMIZE (str):    Optimization algorithm for the inverse problem
    # PREPROCESS (str):  Preprocessing schema for waveform data
    # POSTPROCESS (str): Postprocessing schema for kernels and gradients
    # ==============================================================================
    WORKFLOW: inversion
    SOLVER: specfem3d
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
    #   Verbosity output of SF3 logger. Available from least to most verbosity:
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
    TITLE: docs
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


.. code:: ipython3

    ! tail --lines=54 parameters.yaml


.. parsed-literal::

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
        SCRATCH: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/scratch
        OUTPUT: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/output
        SYSTEM: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/scratch/system
        LOCAL:
        LOGFILE: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/output_sf3.txt
        SOLVER: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/scratch/solver
        SPECFEM_BIN: !!! REQUIRED PATH !!!
        SPECFEM_DATA: !!! REQUIRED PATH !!!
        DATA:
        MASK:
        OPTIMIZE: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/scratch/optimize
        MODEL_INIT: !!! REQUIRED PATH !!!
        MODEL_TRUE:
        FUNC: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/scratch/scratch
        GRAD: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/scratch/evalgrad
        HESS: /home/bchow/REPOSITORIES/seisflows/seisflows/docs/scratch/evalhess


How do I know what parameters need to be set?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   **NOTE**: Required parameters that can not be set to default values
   will be listed as ``!!! REQUIRED PARAMETER !!!``

We can check the required paths and parameters manually by scrolling
through the parameter file, or we can use the
``seisflows par --required`` command to list them out all at once.

.. code:: ipython3

    ! seisflows par --required


.. parsed-literal::

    !!! REQUIRED PARAMETER !!!
    ==========================
    	MATERIALS
    	DENSITY
    	ATTENUATION
    	NT
    	DT
    	FORMAT
    	CASE
    	END
    !!! REQUIRED PATH !!!
    =====================
    	SPECFEM_BIN
    	SPECFEM_DATA
    	MODEL_INIT


Checking parameter validity
---------------------------

You might be asking, how do I know if my parameters are set correctly?
SeisFlows modules feature check() functions which dictate correct
parameter values. You can run ``seisflows init`` to run these check()
functions. Because we have required parameters still left unset in our
parameter file, we expect the ``seisflows init`` function to throw an
error.

.. code:: ipython3

    ! seisflows init


.. parsed-literal::

    ================================================================================
                               PARAMETER FILE READ ERROR                            
                               /////////////////////////                            
    Please check that your parameter file is properly formatted in the YAML format.
    If you have just run 'seisflows configure', you may have some required
    parameters that will need to be filled out before you can proceed. The error
    message is:
    
    could not determine a constructor for the tag 'tag:yaml.org,2002:!'
      in "parameters.yaml", line 147, column 12
    ================================================================================


Let’s set some random variables for the required parameters with the
``seisflows par`` command and try again.

.. code:: ipython3

    ! seisflows par materials elastic
    ! seisflows par density constant
    ! seisflows par attenuation False
    ! seisflows par nt 100
    ! seisflows par dt .05
    ! seisflows par format ascii
    ! seisflows par case data
    ! seisflows par end 1
    ! seisflows par specfem_bin ./
    ! seisflows par specfem_data ./
    ! seisflows par model_init ./


.. parsed-literal::

    MATERIALS: !!! REQUIRED PARAMETER !!! -> elastic
    DENSITY: !!! REQUIRED PARAMETER !!! -> constant
    ATTENUATION: !!! REQUIRED PARAMETER !!! -> False
    NT: !!! REQUIRED PARAMETER !!! -> 100
    DT: !!! REQUIRED PARAMETER !!! -> .05
    FORMAT: !!! REQUIRED PARAMETER !!! -> ascii
    CASE: !!! REQUIRED PARAMETER !!! -> data
    END: !!! REQUIRED PARAMETER !!! -> 1
    SPECFEM_BIN: !!! REQUIRED PATH !!! -> ./
    SPECFEM_DATA: !!! REQUIRED PATH !!! -> ./
    MODEL_INIT: !!! REQUIRED PATH !!! -> ./


.. code:: ipython3

    ! seisflows init


.. parsed-literal::

    instantiating SeisFlows working state in directory: output


Of course we knew that the above parameters were acceptable. But what if
we input an unacceptable parameter into the parameter file and try
again?

.. code:: ipython3

    ! rm -r output/
    ! seisflows par materials visibily_incorrect_value
    ! seisflows init


.. parsed-literal::

    MATERIALS: elastic -> visibily_incorrect_value
    ================================================================================
                                   MODULE CHECK ERROR                               
                                   //////////////////                               
    seisflows.config module check failed with:
    
    solver: MATERIALS must be in ['ELASTIC', 'ACOUSTIC', 'ISOTROPIC', 'ANISOTROPIC']
    ================================================================================


And voila, the module check has thrown an error, and told us (the User)
how to properly set the value of the materials parameter. Hopefully a
combination of thorough explanations in the parameter file section
headers, and error catching with ``seisflows init`` makes crafting your
own parameter file a smooth process.

.. code:: ipython3

    ! head -155 parameters.yaml | tail --lines=38


.. parsed-literal::

    
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
    # FORMAT (float):
    #   Format of synthetic waveforms used during workflow, available options:
    #   ['ascii', 'su']
    # SOURCE_PREFIX (str):
    #   Prefix of SOURCE files in path SPECFEM_DATA. Available ['CMTSOLUTION',
    #   FORCESOLUTION']
    # =============================================================================
    MATERIALS: visibily_incorrect_value
    DENSITY: constant
    ATTENUATION: False
    COMPONENTS: ZNE
    SOLVERIO: fortran_binary
    NT: 100
    DT: .05
    FORMAT: ascii
    SOURCE_PREFIX: CMTSOLUTION


.. code:: ipython3

    ! rm parameters.yaml  # to delete the created file from this working directory
