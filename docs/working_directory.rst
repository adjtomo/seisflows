Working Directory
=================

SeisFlows sets it's own working directory when executing a workflow. Below we explore the working directory set up by the `SPECFEM2D-workstation example <specfem2d_example.html>`__. Working directories may change slightly depending on the chosen workflow, but will more or less follow the same structure.

   **NOTE**: The two SPECFEM2D directories listed below (specfem2d/ &
   specfem2d_workdir/) are not part of a standard SeisFlows working
   directory.

.. code:: ipython3

    %cd /home/bchow/Work/scratch
    ! ls


.. parsed-literal::

    /home/bchow/Work/scratch
    logs	parameters.yaml  sflog.txt    specfem2d
    output	scratch		 sfstate.txt  specfem2d_workdir


--------------

scratch/
--------

The active working directory of SeisFlows where all of the heavy lifting
takes place. Each module in the SeisFlows package may have it’s own
sub-directory where it stores temporary work data. Additionally, we have
two eval*/ directories where objective function evaluation (eval_func)
and gradient evaluation (eval_grad) files are stored.

.. code:: ipython3

    ! ls scratch


.. parsed-literal::

    eval_func  eval_grad  optimize	preprocess  solver  system


.. warning:: 
    As suggested by the name, the scratch/ directory is not for permanent storage, and any data contained within is liable to be changed or removed throughout the workflow. Permanent data storage takes place in the **output/** directory.

solver/
~~~~~~~

A collection of event-specific directories (one directory per event),
each of which is a self contained SPECFEM run directory (i.e., they
contain all the necessary files to run SPECFEM binaries within).

.. note::
    The first (alphabetically) solver in this directory is symlinked as the **mainsolver**. The mainsolver is used to run single-process functions (e.g., gradient smoothing). The mainsolver is also useful for scripting, as the name of the first event may be different from workflow to workflow, so **mainsolver** provides a consistent entry point into the solver subdirectories.

.. code:: ipython3

    ! ls scratch/solver


.. parsed-literal::

    001  002  003  mainsolver


.. code:: ipython3

    ! ls scratch/solver/mainsolver


.. parsed-literal::

    adj_solver.log	combine_vs.log	fwd_solver.log	SEM
    bin		DATA		kernel_paths	traces
    combine_vp.log	fwd_mesher.log	OUTPUT_FILES


.. warning::
    Each solver directory is constructed by copying the `PATH.SPECFEM_BIN` and `PATH.SPECFEM_DATA` directories into each sub-directory. The user must ensure that these directories do not contain large files (e.g., waveform data, large tomography files), otherwise these will be copied N times, where N is the number of events in your workflow. This can quickly run up against storage issues.

The bin/, DATA/ and OUTPUT_FILES/ directories are the same as those
found in SPECFEM. The SEM file defines the locations of the adjoint
sources, which is dictated by SPECFEM. The traces/ directory contains
all of the output waveforms required by this event. They are separated
into observed (obs), synthetic (syn) and adjoint (adj) waveforms.

.. code:: ipython3

    ! ls scratch/solver/mainsolver/traces


.. parsed-literal::

    adj  obs  syn


.. code:: ipython3

    ! ls scratch/solver/mainsolver/traces/obs


.. parsed-literal::

    AA.S0001.BXY.semd


.. code:: ipython3

    # These waveforms are saved into a two-column ASCII format
    ! tail scratch/solver/mainsolver/traces/obs/AA.S0001.BXY.semd


.. parsed-literal::

       251.39999999999998         -1.1814422395268879E-005
       251.45999999999998         -1.1800275583562581E-005
       251.51999999999998         -1.1769315129746346E-005
       251.57999999999998         -1.1721248953632887E-005
       251.63999999999999         -1.1655830825336088E-005
       251.69999999999999         -1.1572872866742356E-005
       251.75999999999999         -1.1472248505521453E-005
       251.81999999999999         -1.1353902449899163E-005
       251.88000000000000         -1.1217847351013855E-005
       251.94000000000000         -1.1064166223014224E-005


optimize/
~~~~~~~~~

Values relating to the optimization algorithm. These variables define
model vectors, misfits, gradient directions and search directions.
Optimization vectors are stored as NumPy arrays and tagged with the .npy
suffix. Optimization scalars are stored as text files and tagged with
the .txt suffix.

Optimization Variable Names are described as:

* m_new: current model vector
* m_old: previous model vector 
* m_try: line search model vector 
* f_new: current objective function value  
* f_old: previous objective function value  
* f_try: line search function value  
* g_new: current gradient direction vector 
* g_old: previous gradient direction vector 
* p_new: current search direction vector 
* p_old: previous search direction vector  

.. code:: ipython3

    ! ls scratch/optimize


.. parsed-literal::

    alpha.txt	f_new.txt  f_try.txt  m_new.npz  output_optim.txt
    checkpoint.npz	f_old.txt  g_old.npz  m_old.npz  p_old.npz


.. code:: ipython3

    import numpy as np
    m_new = np.load("scratch/optimize/m_new.npz")
    print(m_new["vs"])


.. parsed-literal::

    [[3500.0027437  3499.99441921 3499.90777902 ... 3499.77655378
      3499.9021825  3499.99078301]]


.. code:: ipython3

    ! cat scratch/optimize/f_new.txt


.. parsed-literal::

    8.645199999999999153e-04


The ‘checkpoint.npz’ file contains information about the state of the
line search (controlled by the Optimization module). It is used to
resume failed or stopped line searches with minimal redundant use of
computational resources.

.. code:: ipython3

    line_search = np.load("scratch/optimize/checkpoint.npz")
    
    print(vars(line_search)["files"])
    
    print("step count: ", line_search["step_count"])
    print("step lengths: ", line_search["step_lens"])
    print("misfit: ", line_search["func_vals"])


.. parsed-literal::

    ['restarted', 'func_vals', 'step_lens', 'gtg', 'gtp', 'step_count']
    step count:  0
    step lengths:  [0.00000000e+00 2.32268310e+09 3.75818023e+09 1.59087505e+09
     2.82031810e+09]
    misfit:  [0.00127902 0.00086452 0.00172904 0.00259356 0.00345808]


eval_func/ & eval_grad/
~~~~~~~~~~~~~~~~~~~~~~~

Scratch directories containing objective function evaluation and
gradient evaluation files. These include (1) the current **model** being
used for misfit evaluation, and (2) a **residual** file which defines
the misfit for each event. **eval_grad/** also contains **kernels**
which define per-event kernels which are summed and manipulated with the
postprocess module.

.. code:: ipython3

    ! ls scratch/eval_func
    ! echo
    ! ls scratch/eval_grad


.. parsed-literal::

    model  residuals.txt
    
    gradient  kernels  misfit_kernel  model  residuals.txt


.. code:: ipython3

    ! cat scratch/eval_grad/residuals.txt


.. parsed-literal::

    2.41E-02
    2.14E-02
    1.55E-02


.. code:: ipython3

    ! ls scratch/eval_grad/kernels


.. parsed-literal::

    001  002  003


.. code:: ipython3

    ! ls scratch/eval_grad/kernels/001


.. parsed-literal::

    proc000000_bulk_beta_kernel.bin  proc000000_rhop_kernel.bin
    proc000000_bulk_c_kernel.bin	 proc000000_vp_kernel.bin
    proc000000_kappa_kernel.bin	 proc000000_vs_kernel.bin
    proc000000_mu_kernel.bin	 proc000000_weights_kernel.bin
    proc000000_rho_kernel.bin


system & preprocess
~~~~~~~~~~~~~~~~~~~

These two directories are empty in our example problem, but are
catch-all directories where module-specific files can be output. If you
are extending SeisFlows with other base or subclasses, it is preferable
to adhere to this structure where each module only interacts with it’s
own directory.

When ``Pyaflowa`` is chosen as the preprocess module, it stores figures,
log files, and data (in ASDFDataSets) within its scratch directory. It
also specifies parameters for exporting these scratch files to disk for
more permanent storage.

--------------

output/
-------

Output files to be permanently saved (e.g., models, graidents, traces)
can be located in this directory. These are tagged in ascending order.
Because we did not run the finalization task in our SPECFEM2D problem,
the output directory only contains our initial model.

.. code:: ipython3

    ! ls output


.. parsed-literal::

    MODEL_INIT


.. code:: ipython3

    ! ls output/MODEL_INIT


.. parsed-literal::

    proc000000_vp.bin  proc000000_vs.bin


--------------

logs/
-----

Where any text logs are stored. If running on a cluster, all submitted
jobs will be instructed to write their logs into this directory.
Additionally, if a workflow is resumed (previous log files exist in the
other directory) copies are saved to this directory.

.. code:: ipython3

    ! ls logs


.. parsed-literal::

    0001_00.log  0002_02.log  0004_01.log  0006_00.log	    parameters_002.yaml
    0001_01.log  0003_00.log  0004_02.log  0006_01.log	    parameters_003.yaml
    0001_02.log  0003_01.log  0005_00.log  0006_02.log	    sflog_001.txt
    0002_00.log  0003_02.log  0005_01.log  0007_00.log	    sflog_002.txt
    0002_01.log  0004_00.log  0005_02.log  parameters_001.yaml  sflog_003.txt


--------------

sflog.txt
---------

The main log file for SeisFlows, where all log statements written to
stdout are recorded during a workflow. Allows a user to come back to a
workflow and understand the tasks completed and any important
information collected during the workflow

.. code:: ipython3

    ! head -50 sflog.txt


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


--------------

sfstate.txt
-----------

A state file which tracks the progress of a workflow, allowing the User
to quickly resumed stopped or failed workflows without wasting
computational resources. The State file simply contains the names of
functions contained in the Workflow task list, as well as their
respective status, which can be ‘completed’, ‘failed’, or not available.

.. code:: ipython3

    ! cat sfstate.txt


.. parsed-literal::

    # SeisFlows State File
    # Tue Aug 16 14:33:17 2022
    # Acceptable states: 'completed', 'failed'
    # =======================================
    evaluate_initial_misfit: completed
    run_adjoint_simulations: completed
    postprocess_event_kernels: completed
    evaluate_gradient_from_kernels: completed
    initialize_line_search: completed
    perform_line_search: completed
    iteration: 1

When submitting a workflow with an existing state file, the workflow
will check the status of each function. ‘Completed’ functions will be
skipped over. ‘Failed’ functions will be re-run. Users can delete lines
from the state file or change status’ manually to re-run tasks within
the list, taking care about the current configuration of the working
directory, which is intrinsically tied to the task list.

For ‘Inversion’ workflows, the current ‘Iteration’ is also saved,
meaning re-submitted workflows will start at the previously checkpointed
iteration.
