Working Directory Structure
===========================

SeisFlows3 hardcodes it’s own working directory when executing a
workflow. Below we explore the working directory set up by the
SPECFEM2D-workstation example. Working directories may change slightly
depending on the chosen workflow, but will more or less follow the
following structure. The two specfem2d directories listed below are not
part of the SeisFlows3 working directory.

.. code:: ipython3

    %cd ~/Work/official/workshop_pyatoa_sf3/ex1_specfem2d_workstation
    ! ls


.. parsed-literal::

    /home/bchow/Work/official/workshop_pyatoa_sf3/ex1_specfem2d_workstation
    logs	output_sf3.txt	 scratch	    stats
    output	parameters.yaml  specfem2d_workdir


--------------

scratch/
--------

The active working directory of SeisFlows3 where all of the heavy
lifting takes place. Each module in the SeisFlows3 package may have it’s
own sub-directory where it stores temporary work data. Additionally, we
have two eval*/ directories where objective function evaluation
(evalfunc) and gradient evaluation (evalgrad) files are stored.

.. code:: ipython3

    ! ls scratch


.. parsed-literal::

    evalfunc  evalgrad  optimize  preprocess  solver  system


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

    bin  DATA  kernel_paths  mesher.log  OUTPUT_FILES  SEM	solver.log  traces


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

    alpha.npy  f_old.txt  g_old.npy  m_new.npy  p_old.npy
    f_new.txt  f_try.txt  LBFGS	 m_old.npy


.. code:: ipython3

    import numpy as np
    m_new = np.load("scratch/optimize/m_new.npy")
    print(m_new)


.. parsed-literal::

    [5800.         5800.         5800.         ... 3499.77655379 3499.9021825
     3499.99078301]


.. code:: ipython3

    ! cat scratch/optimize/f_new.txt


.. parsed-literal::

    2.591424e-03


evalfunc/ & evalgrad/
~~~~~~~~~~~~~~~~~~~~~

Scratch directories containing objective function evaluation and
gradient evaluation files. These include (1) the current **model** being
used for misfit evaluation, and (2) **residuals** which define the
misfit for each event. **evalgrad/** also contains **kernels** which
define per-event kernels which are summed and manipulated with the
postprocess module.

.. code:: ipython3

    ! ls scratch/evalfunc
    ! echo
    ! ls scratch/evalgrad


.. parsed-literal::

    model  residuals
    
    kernels  model	residuals


.. code:: ipython3

    ! ls scratch/evalgrad/residuals


.. parsed-literal::

    001  002  003


.. code:: ipython3

    ! cat scratch/evalgrad/residuals/001


.. parsed-literal::

    2.413801941841247842e-02
    2.413801941841247842e-02
    2.413801941841247842e-02


.. code:: ipython3

    ! ls scratch/evalgrad/kernels


.. parsed-literal::

    001  002  003  sum


.. code:: ipython3

    ! ls scratch/evalgrad/kernels/sum


.. parsed-literal::

    proc000000_vp_kernel.bin  proc000000_vs_kernel.bin


system & preprocess
~~~~~~~~~~~~~~~~~~~

These two directories are empty in our example problem, but are
catch-all directories where module-specific files can be output. If you
are extending SeisFlows3 with other base or subclasses, it is preferable
to adhere to this structure where each module only interacts with it’s
own directory

.. code:: ipython3

    ! ls scratch/system
    ! ls scratch/preprocess

--------------

output/
-------

The current active state of SeisFlows3, containing pickle (.p) and JSON
files which describe a Python environment of a current workflow.
Additionally files to be permanently saved (e.g., models, graidents,
traces) can be located here. These are tagged in ascending order, e.g.,
model_0001 refers to the updated model derived during the first
iteration.

.. code:: ipython3

    ! ls output


.. parsed-literal::

    gradient_0001  seisflows_optimize.p	  seisflows_solver.p
    kwargs	       seisflows_parameters.json  seisflows_system.p
    model_0001     seisflows_paths.json	  seisflows_workflow.p
    model_init     seisflows_postprocess.p
    model_true     seisflows_preprocess.p


.. code:: ipython3

    ! ls output/model_0001


.. parsed-literal::

    proc000000_vp.bin  proc000000_vs.bin


.. code:: ipython3

    ! ls output/gradient_0001


.. parsed-literal::

    proc000000_vp_kernel.bin  proc000000_vs_kernel.bin


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

    output_sf3_001.txt  parameters_001.yaml


--------------

stats/
------

Text files describing the optimization statistics of the current
workflow. This directory is only relevant if you are running an
inversion workflow.

.. code:: ipython3

    ! ls stats


.. parsed-literal::

    factor.txt	      line_search.txt  slope.txt	theta.txt
    gradient_norm_L1.txt  misfit.txt       step_count.txt
    gradient_norm_L2.txt  restarted.txt    step_length.txt


.. code:: ipython3

    ! cat stats/step_count.txt


.. parsed-literal::

    ITER          STEP_COUNT
    ====  ==================
       1        0.000000E+00
       1        2.000000E+00


--------------

output_sf3.txt
--------------

The main log file for SeisFlows3, where all log statements written to
stdout are recorded during a workflow.

.. code:: ipython3

    ! head -50 output_sf3.txt


.. parsed-literal::

    2022-04-29 16:45:35 | initializing SeisFlows3 in sys.modules
    2022-04-29 16:45:39 | copying par/log file to: /home/bchow/Work/official/workshop_pyatoa_sf3/ex1_specfem2d_workstation/logs/output_sf3_001.txt
    2022-04-29 16:45:39 | copying par/log file to: /home/bchow/Work/official/workshop_pyatoa_sf3/ex1_specfem2d_workstation/logs/parameters_001.yaml
    2022-04-29 16:45:39 | exporting current working environment to disk
    2022-04-29 16:45:39 | 
    ////////////////////////////////////////////////////////////////////////////////
                       WORKFLOW WILL STOP AFTER FUNC: 'finalize'                    
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 16:45:39 | 
    ================================================================================
                              STARTING INVERSION WORKFLOW                           
    ================================================================================
    2022-04-29 16:45:39 | 
    ////////////////////////////////////////////////////////////////////////////////
                                    ITERATION 1 / 1                                 
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 16:45:39 | 
    ////////////////////////////////////////////////////////////////////////////////
                                PERFORMING MODULE SETUP                             
    ////////////////////////////////////////////////////////////////////////////////
    2022-04-29 16:45:39 | misfit function is: 'waveform'
    2022-04-29 16:45:40 | writing line search history file:
    /home/bchow/Work/official/workshop_pyatoa_sf3/ex1_specfem2d_workstation/stats/line_search.txt
    2022-04-29 16:45:40 | checking poissons ratio for: 'm_new.npy'
    2022-04-29 16:45:40 | model parameters (m_new.npy i01s00):
    2022-04-29 16:45:40 | 5800.00 <= vp <= 5800.00
    2022-04-29 16:45:40 | 3500.00 <= vs <= 3500.00
    2022-04-29 16:45:40 | 0.21 <= pr <= 0.21
    2022-04-29 16:45:41 | setting up solver on system...
    2022-04-29 16:45:41 | checkpointing working environment to disk
    2022-04-29 16:45:42 | exporting current working environment to disk
    2022-04-29 16:45:43 | running task solver.setup 3 times
    2022-04-29 16:45:43 | initializing 3 solver directories
    2022-04-29 16:45:50 | source 001 symlinked as mainsolver
    2022-04-29 16:45:50 | generating 'data' with MODEL_TRUE synthetics
    2022-04-29 16:45:57 | running mesh generation for MODEL_INIT
    2022-04-29 16:46:27 | 
    ================================================================================
                                 INITIALIZING INVERSION                             
    ================================================================================
    2022-04-29 16:46:27 | 
    EVALUATE OBJECTIVE FUNCTION
    --------------------------------------------------------------------------------
    2022-04-29 16:46:27 | saving model 'm_new.npy' to:
    /home/bchow/Work/official/workshop_pyatoa_sf3/ex1_specfem2d_workstation/scratch/evalgrad/model
    2022-04-29 16:46:28 | evaluating objective function 3 times on system...
    2022-04-29 16:46:28 | checkpointing working environment to disk
    2022-04-29 16:46:29 | exporting current working environment to disk
    2022-04-29 16:46:30 | running task solver.eval_func 3 times
    2022-04-29 16:46:30 | running forward simulations

