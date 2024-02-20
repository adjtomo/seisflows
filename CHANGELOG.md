### Major Changes
- Workflow:
    - State file changed states from words to integers (completed -> 1, 
      failed -> -1, pending -> 0), and full state file is created at workflow
      setup, rather than one by one as each function completes
    - Inversion:
        - **Changed total misfit summation from summed + squared to L1 norm.**
        - Line search has been broken into multiple functions to facilitate
          restarting failed line searches: 
          perform_line_search() -> evaluate_line_search_misfit() + update_line_search()
    - NoiseInversion:
        - **Modifies the Inversion class to invert for ZZ, RR and TT empirical 
          Green's functions**
        - Utility function to rotate EE, EN, NE, NN waveforms to RR and TT
        - Utility function to rotate RR and TT adjoint sources to EE, EN, NE, NN
        - Utility function to convert STATIONS file to N SOURCE files since 
          virtual sources are required for ambient noise adjoint tomography
        - Functionality to create a source-receiver lookup table that contains
          information on azimuth + backazimuth required for RR/TT kernels
- Preprocessing:
    - Default (major upgrades)
        - Force 'obs' data to match 'syn' data sampling rate
        - **Parallelized misfit quantification with concurrent futures**
        - Zero'd adjoint source generation now occurs at module setup 
          (parallelized)
        - Added additional normalization options
        - Allows selection by component for misfit quantification
        - Improved obs-syn file match validation
        - Removed the ability to `sum_residuals`, which required Preprocess to 
          know too many things about the workflow. Now handled by Workflow.
        - **Added a simple waveform plotter to show obs, syn and adjoint source**
        - Split off preprocessing functions (mute, normalize) to `tools` which
          Preprocess can import
    - Pyaflowa
        - Removed `client` parameter to match Pyatoa > 0.3.0
        - Allow processing only for specific components
        - Data reading abstraction simplified, no longer builds paths from parts
          but instead explcitely reads data + metadata like Default preproc.
        - **Pyflex preset now directly part of parameter file so that User can
          edit them directly**
- System:
    - General:
        - `tasktime` now set in the top parent class
        - Allow custom tasktimes for functions run with `System.run`
        - **`rerun` function tells System to re-run failed jobs some number of 
          times to deal with randomly failing tasks that usually work once you
          run them again**; added to TestFlow
    - Cluster (and derived classes):
        - **New parameter `array`: For debug purposes, allow running only specific
          task IDs to e.g., re-run failed processes. Input style follows SLURM
          array argument** 
        - **Submit jobs directly to the login node with the -l/--login flag**
        - Non-zero exit code error catching added to concurrent future calls
        - **Overhauled job monitoring system. Notably, does not break on first  
          job failure, but rather to wait until all jobs are finished. Tied into System `rerun` feature**
    - Slurm (and derived classes): 
        - Added a timeout counter and extended timeout value for checking
          output of `sacct` for queue checking due to premature job exits with
          empty `sacct` returns (i.e., it takes a while for compute nodes to 
          spin up and be visible in `sacct`)
- Optimize:
    - Major reorganization, breaking major monolithic functions into smaller pieces
    - Improved function and checkpoint order with the intent of easier failure recovery
    - Internal restart and revert functions for manually stepping back failed line searches
- Solver: 
    - API change: solver.combine() made more generic and no longer hardcodes
      assumed directory structure
    - Parameter change: `density` -> `update_density`
    - Model parameter checks removed from Solver's abilities. These are now 
      handled by Workflow
    - Takes over responsibility for renaming adjoint sources 
    - Takes over responsibility for obs-syn filename matching prior to preproc.
- Command Line Tool:
    - `seisflows submit --login` -> `seisflows submit --direct` for submitting
      your workflow directly to the login/home node and not to a cluster node
    - `seisflows configure` makes more clear which paths are default/not 
      important, and which paths are required
    - `seisflows setup` <-> `seisflows init` namespace change as the names make
      more sense in this order. `init` starts a blank working directory, `setup`
      runs module setup functions (like directory creation)
- New Dependencies: 
    - PyPDF: for PDF mergers in Pyaflowa preprocessing
    - PySEP: for SPECFEM-specific read functions 

### Minor Changes
- Solver: 
    - Parallelized directory initialization w/ concurrent futures
    - Kernel renaming defined as a separate function (previously part of 
      adjoint simulation), so that it can be called by debugger
- Optimization:
    - Improves step length overwrite log messaging
    - Skips initial line search calculation if initial step length requested
- Workflow: Allow other workflows to overwrite the default location where 
  synthetic waveforms are saved
- System: improved file system management by organizing spawned process log files
  and removing scratch/system directory each iteration
- Preprocessing Pyaflowa: improved output file management, export and removal 
- Model checking now occurs in Workflow rather than Solver functions
- **Workflow data preparation now symlinks in real data rather than copies it,
  to avoid heavy file overhead**
- **Removed unnused `graphics` utilities and adopted all of Pyatoa's PNG and PDF
  image manipulation utilities for use in Pyaflowa preprocessing (adds PyPDF as
  dependency)**
- Removed large sections of commented out code from command line tool
- Single source version number in `pyproject.toml`
- **Logger aesthetic change to show first four letters of message type rather than
  first letter (e.g., I -> INFO, W -> WARN, D -> DEBU)**
- **Model parameter check now includes mean values in addition to min and max**
- **New CLI tool: `seisflows print tasks` to get source names and relevant task ID**
- Model tool now breaks on check if any NaNs present in model arrays
- `Optimize`: split up some internal functions for easier separation of tasks 
  that were previously all mashed together

### Bugfixes
- **Major: SPECFEM3D_GLOBE based solvers were NOT updating the model during 
  Inversion workflows ecause `xmeshfem3D` was not being called, and therefore 
  not updating database files.**
- **Pyaflowa was checking the incorrect values for windows returned from Manager 
  causing incorrect misfit calculation and an inability for the inversion to reduce misfit**
- **Cross-correlation Traveltime misfit function was not squared, allowing CC 
  values to be negative. Now follows Tromp (2005) where we square the time shift**
- `mpiexec` was being set inside System initiation, causing check statement to 
  fail quietly
- System.Cluster.Run now passes User-defined arguments for log level and 
  verbosity to each child process allowing for uniform logs for all jobs
- `seisflows swap` allow paths to be set relative or absolute, previously they
  were forced to absolute
- Old Optimization files (e.g., m\_old) were not being deleted due to missing 
  file extensions. Not critical because they were not used, and overwritten
- **`seisflows setpar` was not properly setting FORTRAN double precision values. 
  Added some better catches for `setpar` as it was quietly failing when files
  were nonexistent**
- LBFGS line search restart `step_count_max` was not being evaluated properly
- `seisflows configure` will no longer try to configure an already configured   
  file
- **Concurrency: better all-around error catching for any functions that are 
  parallelized by concurrent futures. Previously these functions failed quietly**
- solver.specfem3d_globe was not recognizing custom model types

### Misc.
- Removed hard requirement that `import_seisflows` required all Workflows have
  `modules` as their first argument. Only Forward workflow requires.
- Removed Optimize load checkpoint from inversion setup because it was already
  run by Optimize setup
