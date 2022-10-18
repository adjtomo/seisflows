Parameter File
==============

The parameter file is the central control object for a SeisFlows
workflow. Here we take a look at the anatomy of a parameter file.
Parameter files in SeisFlows are formatted in the `YAML format (YAML
Ain’t Markup Language) <https://pyyaml.org/wiki/PyYAMLDocumentation>`__.

Template
--------

Each workflow starts with the module-only template parameter file which
defines the core modules of the package. Your choices for each of these
modules will determine which paths and parameters are included in the
full parameter file. Running ``seisflows setup`` from the command line
will create the template file.

.. code:: ipython3

    ! seisflows setup -h


.. parsed-literal::

    usage: seisflows setup [-h] [-f]
    
    In the specified working directory, copy template parameter file containing
    only module choices, and symlink source code for both the base and super
    repositories for easy edit access. If a parameter file matching the provided
    name exists in the working directory, a prompt will appear asking the user if
    they want to overwrite.
    
    optional arguments:
      -h, --help   show this help message and exit
      -f, --force  automatically overwrites existing parameter file


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


How do I choose modules?
~~~~~~~~~~~~~~~~~~~~~~~~

As seen above, each of the modules comes with a default value which
represents the base class\* for this module.

* For an explanation of base classes and Python inheritance, see the `inheritance page <inheritance.html>`__ 

These default values are likely not suitable for all, e.g., if you want
to run an inversion and not a forward workflow, or use SPECFEM3D not
SPECFEM2D. To see all available module options, use the
``seisflows print modules`` command.

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


How do I change modules?
~~~~~~~~~~~~~~~~~~~~~~~~

Feel free to use any text editor, or use the ``seisflows par`` command
to make changes directly from the command line. For example, say we want
to use SPECFEM3D as our solver module.

This is also covered in the `command line tool page <command_line_tool.html>`__

.. code:: ipython3

    # Changes the current parameter to the given value
    ! seisflows par solver specfem3d


.. parsed-literal::

    solver: specfem2d -> specfem3d


.. code:: ipython3

    # Prints out the current parameter value
    ! seisflows par solver


.. parsed-literal::

    solver: specfem3d


How do I create a full parameter file?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module-only parameter file serves as as a template for dynamically
generating the full parameter file. Since each module requires it’s own
unique set of parameters and paths, each parameter file will look
different. We use the ``seisflows configure`` command to complete the
file.

.. code:: ipython3

    ! seisflows configure -h


.. parsed-literal::

    usage: seisflows configure [-h] [-a]
    
    SeisFlows parameter files will vary depending on chosen modules and their
    respective required parameters. This function will dynamically traverse the
    source code and generate a template parameter file based on module choices.
    The resulting file incldues docstrings and type hints for each parameter.
    Optional parameters will be set with default values and required parameters
    and paths will be marked appropriately. Required parameters must be set before
    a workflow can be submitted.
    
    optional arguments:
      -h, --help            show this help message and exit
      -a, --absolute_paths  Set default paths relative to cwd


.. code:: ipython3

    ! seisflows configure

Below we will take a look at the parameter file we just created

Anatomy of a parameter file
---------------------------

Each of SeisFlows’ modules will define its own section in the parameter
file, separated by a header of comments representing the docstring.
Within each header, parameter names, types and descriptions are listed.
At the bottom of the parameter file, there is a section defining paths
required by SeisFlows. Section headers will look something:

.. code:: ipython3

    # =============================================================================
    # MODULE
    # ------
    # Module description 
    #
    # Parameters
    # ----------
    # :type parameter: type
    # :param paramter: description
    # ...
    # =============================================================================
    parameter: value

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
    # workflow (str):    The types and order of functions for running SeisFlows
    # system (str):      Computer architecture of the system being used
    # solver (str):      External numerical solver to use for waveform simulations
    # preprocess (str):  Preprocessing schema for waveform data
    # optimize (str):    Optimization algorithm for the inverse problem
    # ==============================================================================
    workflow: forward
    system: workstation
    solver: specfem3d
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
    #        not attempt to generate data.
    #    :type export_traces: bool
    #    :param export_traces: export all waveforms that are generated by the
    #        external solver to `path_output`. If False, solver traces stored in
    #        scratch may be discarded at any time in the workflow
    #    :type export_residuals: bool
    #    :param export_residuals: export all residuals (data-synthetic misfit) that
    #        are generated by the external solver to `path_output`. If False,
    #        residuals stored in scratch may be discarded at any time in the workflow
    #
    #        
    # =============================================================================
    data_case: data
    export_traces: False
    export_residuals: False
    # =============================================================================
    #
    #    Workstation System
    #    ------------------
    #    Runs tasks in serial on a local machine.
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


.. code:: ipython3

    ! tail parameters.yaml


.. parsed-literal::

    path_model_true: null
    path_state_file: /Users/Chow/Repositories/seisflows/docs/notebooks/sfstate.txt
    path_data: null
    path_par_file: /Users/Chow/Repositories/seisflows/docs/notebooks/parameters.yaml
    path_log_files: /Users/Chow/Repositories/seisflows/docs/notebooks/logs
    path_output_log: /Users/Chow/Repositories/seisflows/docs/notebooks/sflog.txt
    path_specfem_bin: null
    path_specfem_data: null
    path_solver: /Users/Chow/Repositories/seisflows/docs/notebooks/scratch/solver
    path_preconditioner: null


How do I know how parameters need to be set?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most SeisFlows parameters come with reasonable default values. The
docstrings headers will also list the expected type and available
options (if any). You may also run the ``seisflows check`` command which
verifies that parameters are set correctly.

.. code:: ipython3

    ! seisflows check


.. parsed-literal::

    
    ================================================================================
                                    PARAMETER ERRROR                                
                                    ////////////////                                
    `path_specfem_bin` must exist and must point to directory containing SPECFEM
    executables
    ================================================================================


.. code:: ipython3

    ! rm parameters.yaml  # to delete the created file from this working directory
