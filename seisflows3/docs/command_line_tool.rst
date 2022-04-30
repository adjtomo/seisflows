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


Setting up a working directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    ! seisflows setup


.. parsed-literal::

    creating parameter file: parameters.yaml


.. code:: ipython3

    ! ls


.. parsed-literal::

    parameters.yaml


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
    #       > seisflows print module
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


