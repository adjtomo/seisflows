Commmand Line Tool
==================

``SeisFlows3`` is primarily interacted with via command line calls and a
parameter file. In this page we explain how to use this command line
tool to create a SeisFlows3 parameters file, edit and configure it, and
establish a SeisFlows3 working directory. We also provide explanation
for other command line options which act as helper utilities for
improved package control.

After installing SeisFlows3 into a Conda environment, the ``seisflows``
command will be available directly from the command line. To access the
help dialogue, you can type ``seisflows`` or ``seisflows -h``

.. code:: ipython3

    ! seisflows


.. parsed-literal::

    usage: seisflows [-h] [-w [WORKDIR]] [-p [PARAMETER_FILE]]
                     {setup,configure,init,submit,resume,restart,clean,par,sempar,check,print,convert,reset,debug,edit,examples}
                     ...
    
    ================================================================================
    
                         SeisFlows3: Waveform Inversion Package                     
    
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
        par                 View and edit SeisFlows3 parameter file
        sempar              View and edit SPECFEM parameter file
        check               Check state of an active environment
        print               Print information related to an active environment
        convert             Convert model file format
        reset               Reset modules within an active state
        debug               Start interactive debug environment
        edit                Open source code file in text editor
        examples            Look at and run pre-configured example problems
    
    'seisflows [command] -h' for more detailed descriptions of each command.


Setting up a parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

seisflows setup
^^^^^^^^^^^^^^^

The first step of any SeisFlows3 workflow is to set up a working
directory, which begins by establishing a blank parameter file. The
``seisflows setup`` command copies in a template parameter file. Ideally
your working directory will be empty to avoid file conflicts.

.. code:: ipython3

    %cd ~/Work/scratch
    ! ls


.. parsed-literal::

    /home/bchow/Work/scratch


.. code:: ipython3

    ! seisflows setup -h


.. parsed-literal::

    usage: seisflows setup [-h] [-s] [-f]
    
    In the specified working directory, copy template parameter file containing
    only module choices, and symlink source code for both the base and super
    repositories for easy edit access. If a parameter file matching the provided
    name exists in the working directory, a prompt will appear asking the user if
    they want to overwrite.
    
    optional arguments:
      -h, --help     show this help message and exit
      -s, --symlink  symlink source code into the working directory
      -f, --force    automatically overwrites existing parameter file


.. code:: ipython3

    # The '-f' flag (force) will overwrite any existing parameter file
    ! seisflows setup -f


.. parsed-literal::

    creating parameter file: parameters.yaml


Having a look at the template parameters.yaml file that was just
generated, we can see that it contains some pre-defined default values
for the core SeisFlows3 modules. Each of these modules defines it’s own
set of unique parameters which make up a workflow.

.. code:: ipython3

    ! ls
    ! wc -l parameters.yaml  # List the number of lines in the file


.. parsed-literal::

    parameters.yaml
    32 parameters.yaml


.. code:: ipython3

    ! cat parameters.yaml


.. parsed-literal::

    # //////////////////////////////////////////////////////////////////////////////
    #
    #                        SeisFlows3 YAML Parameter File
    #
    # //////////////////////////////////////////////////////////////////////////////
    #
    # Modules correspond to the structure of the source code, and determine
    # SeisFlows3' behavior at runtime. Each module requires its own sub-parameters.
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
    # WORKFLOW (str):    The method for running SeisFlows3; equivalent to main()
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


seisflows configure
^^^^^^^^^^^^^^^^^^^

We can now run the ``seisflows configure`` command which will build out
our parameter file based on the module choices provided in the parameter
file.

.. code:: ipython3

    ! seisflows configure -h


.. parsed-literal::

    usage: seisflows configure [-h] [-r]
    
    SeisFlows parameter files will vary depending on chosen modules and their
    respective required parameters. This function will dynamically traverse the
    source code and generate a template parameter file based on module choices.
    The resulting file incldues docstrings and type hints for each parameter.
    Optional parameters will be set with default values and required parameters
    and paths will be marked appropriately. Required parameters must be set before
    a workflow can be submitted.
    
    optional arguments:
      -h, --help            show this help message and exit
      -r, --relative_paths  Set default paths relative to cwd


.. code:: ipython3

    ! seisflows configure


.. parsed-literal::

    filling parameters.yaml w/ default values


.. code:: ipython3

    ! head -200 parameters.yaml | tail -n 82  # have a look at the middle of the file
    ! echo
    ! wc -l parameters.yaml


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
    #   seisflows3.plugins.preconds for available choices
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
    
    306 parameters.yaml


We can see that our parameter file is over 300 lines now, too cumbersome
to print on the page. The length of the file mostly arises from the
header, as each parameter gets it’s own entry with the parameter’s type,
docstring, and any available options.

Parameters that are required by the workflow but do not come with
pre-set default values will be labelled with
``!!! REQUIRED PARAMETER !!!``. Similarly required path definitions,
which come at the end of the file, are labelled with the
``!!! REQUIRED PATH !!!`` value.

Filling out the parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seisflows par
^^^^^^^^^^^^^

It’s easy enough to open your favorite text editor to make adjustments
to the parameter file, however the ``seisflows par`` command makes
things easier by allowing you to view and edit values from the command
line. This makes it convenient to change parameters, and also allows you
to script your workflow setup for improved reproducibility.

.. code:: ipython3

    ! seisflows par -h


.. parsed-literal::

    usage: seisflows par [-h] [-p] [-r] [parameter] [value]
    
    Directly edit values in the parameter file by providing the parameter and
    corresponding value. If no value is provided, will simply print out the
    current value of the given parameter. Works also with path names.
    
    positional arguments:
      parameter         Parameter to edit or view, (case independent).
      value             Optional value to set parameter to. If not given, will
                        print out current parameter. If given, will replace
                        current parameter with new value. Set as 'null' for
                        NoneType and set '' for empty string
    
    optional arguments:
      -h, --help        show this help message and exit
      -p, --skip_print  Skip the print statement which is typically sent to stdout
                        after changing parameters.
      -r, --required    Only list parameters which have not been set as a default
                        value, typically set with some attention catching
                        argument. 'parameter' and 'value' will be ignored.


The -r (–required) flag tells us which parameters need to be set by the
user

.. code:: ipython3

    ! seisflows par -r


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


We can view (but not modify) parameters by giving a single argument to
the par command

.. code:: ipython3

    ! seisflows par end


.. parsed-literal::

    END: !!! REQUIRED PARAMETER !!!


and we can edit the given parameter by providing a second argument to
the par command

.. code:: ipython3

    ! seisflows par end 1


.. parsed-literal::

    END: !!! REQUIRED PARAMETER !!! -> 1


seisflows sempar
^^^^^^^^^^^^^^^^

The ``seisflows sempar`` command behaves the same as the ``par``
command, except is used to edit a SPECFEM2D/3D/3D_GLOBE Par_file. It has
the same call structure as ``par``.

Setting up an active working state
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An active SeisFlows3 working state is simply a Python environment with
the SeisFlows3 library defined based on the given parameter file. In
order to establish a working state, we need to set all required paths
and parameters. We can look at the parameter file header to determine
valid options for each parameter.

.. code:: ipython3

    ! head -130 parameters.yaml | tail -n 10  


.. parsed-literal::

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


.. code:: ipython3

    # We use the `-p` flag to turn off stdout printing
    ! seisflows par materials elastic -p
    ! seisflows par density constant -p
    ! seisflows par attenuation False -p
    ! seisflows par nt 100 -p
    ! seisflows par dt .01 -p
    ! seisflows par f0 .5 -p
    ! seisflows par format ascii -p
    ! seisflows par case synthetic -p
    
    # Required paths can similarly be set the `par` command
    ! seisflows par specfem_bin ./ -p
    ! seisflows par specfem_data ./ -p
    ! seisflows par model_init ./ -p

seisflows init
^^^^^^^^^^^^^^

To initiate a working state, we run ``seisflows init``. This registers
the parameter file into Python’s sys.modules. It runs parameter check
functions to ensure that parameters have been set correctly, and then
saves the active working state as a set of pickle (.p) files which can
be used to resume active workflows.

.. code:: ipython3

    ! seisflows init


.. parsed-literal::

    
    ================================================================================
                                   MODULE CHECK ERROR                               
                                   //////////////////                               
    seisflows.config module check failed with:
    
    workflow: CASE == SYNTHETIC requires PATH.MODEL_TRUE
    ================================================================================


Oops, as we can see the parameter check has caught that a given
parameter requires a certain path to be set which is currently blank.
Let’s amend and try again

.. code:: ipython3

    ! seisflows par model_true ./ -p
    ! seisflows init


.. parsed-literal::

    instantiating SeisFlows3 working state in directory: output


.. code:: ipython3

    ! ls
    ! echo
    ! ls output


.. parsed-literal::

    output	parameters.yaml
    
    seisflows_optimize.p	   seisflows_postprocess.p  seisflows_system.p
    seisflows_parameters.json  seisflows_preprocess.p   seisflows_workflow.p
    seisflows_paths.json	   seisflows_solver.p

