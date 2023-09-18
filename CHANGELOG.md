# SeisFlows Changelog

## v3.0.0

Ambient Noise Adjoint Tomography (ANAT) implementation and associated package
changes.

### Major Changes
- Workflow:
    - State file changed states from words to integers (completed -> 1, 
      failed -> -1, pending -> 0), and full state file is created at workflow
      setup, rather than one by one as each function completes
    - Inversion:
        - Changed total misfit summation from summed + squared to L1 norm.
        - Line search has been broken into multiple functions to facilitate
          restarting failed line searches: 
          perform_line_search() -> evaluate_line_search_misfit() + update_line_search()
    - NoiseInversion:
        - Modifies the Inversion class to invert for ZZ, RR and TT empirical 
          Green's functions
        - Utility function to rotate EE, EN, NE, NN waveforms to RR and TT
        - Utility function to rotate RR and TT adjoint sources to EE, EN, NE, NN
        - Utility function to convert STATIONS file to N SOURCE files since 
          virtual sources are required for ambient noise adjoint tomography
        - Functionality to create a source-receiver lookup table that contains
          information on azimuth + backazimuth required for RR/TT kernels
- Preprocessing:
    - Default (major upgrades)
        - Force 'obs' data to match 'syn' data sampling rate
        - Parallelized misfit quantification with concurrent futures
        - Zero'd adjoint source generation now occurs at module setup 
          (parallelized)
        - Added additional normalization options
        - Allows selection by component for misfiti quantification
        - Improved obs-syn file match validation
        - Removed the ability to `sum_residuals`, which required Preprocess to 
          know too many things about the workflow. Now handled by Workflow.
        - Added a simple waveform plotter to show obs, syn and adjoint source
        - Split off preprocessing functions (mute, normalize) to `tools` which
          Preprocess can import
    - Pyaflowa
        - Removed `client` parameter to match Pyatoa > 0.3.0
        - Allow processing only for specific components
        - Data reading abstraction simplified, no longer builds paths from parts
          but instead explcitely reads data + metadata like Default preproc.
        - Pyflex preset now directly part of parameter file so that User can
          edit them directly
- System:
    - Cluster (and derived classes):
        - New parameter `array`: For debug purposes, allow running only specific
          task IDs to e.g., re-run failed processes. Input style follows SLURM
          array argument 
        - Submit jobs directly to the login node with the -l/--login flag
        - Non-zero exit code error catching added to concurrent future calls
    - Slurm (and derived classes): 
        - Added a timeout counter and extended timeout value for checking
          output of `sacct` for queue checking due to premature job exits with
          empty `sacct` returns (i.e., it takes a while for compute nodes to 
          spin up and be visible in `sacct`)
- Solver: 
    - API change: solver.combine() made more generic and no longer hardcodes
      assumed directory structure
    - Parameter change: `density` -> `update_density`
    - Model parameter checks removed from Solver's abilities. These are now 
      handled by Workflow
    - Takes over responsibility for renaming adjoint sources 
    - Takes over responsibility for obs-syn filename matching prior to preproc.
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
- Model checking now occurs in Workflow rather than Solver functions
- Workflow data preparation now symlinks in real data rather than copies it,
  to avoid heavy file overhead
- Removed unnused `graphics` utilities and adopted all of Pyatoa's PNG and PDF
  image manipulation utilities for use in Pyaflowa preprocessing (adds PyPDF as
  dependency)
- Removed large sections of commented out code from command line tool
- Single source version number in `pyproject.toml`
- Logger aesthetic change to show first four letters of message type rather than
  first letter (e.g., I -> INFO, W -> WARN, D -> DEBU)
- Model parameter check now includes mean values in addition to min and max

### Bugfixes
- Major: SPECFEM3D\_GLOBE based solvers were NOT updating the model during 
  Inversion workflows ecause `xmeshfem3D` was not being called, and therefore 
  not updating database files.
- Cross-correlation Traveltime misfit function was not squared, allowing CC 
  values to be negative. Now follows Tromp (2005) where we square the time shift
- `mpiexec` was being set inside System initiation, causing check statement to 
  fail quietly
- System.Cluster.Run now passes User-defined arguments for log level and 
  verbosity to each child process allowing for uniform logs for all jobs
- `seisflows swap` allow paths to be set relative or absolute, previously they
  were forced to absolute
- Old Optimization files (e.g., m\_old) were not being deleted due to missing 
  file extensions. Not critical because they were not used, and overwritten

### Misc.
- Removed hard requirement that `import_seisflows` required all Workflows have
  `modules` as their first argument. Only Forward workflow requires.
- Removed Optimize load checkpoint from inversion setup because it was already
  run by Optimize setup



## v2.1.1

- Updates and simplifies install procedure using 'environment.yml' and 
  'pyproject.toml' files. 
- Docs: Adds contributor's guide message to main docs page

## v2.2.0

### #155
- Synthetics and observations can have a different format
- New parameters `obs_data_format` and `syn_data_format` replace old parameter
  `data_format`
- Observations can have SAC format
- Added <unit_output> variable to default preprocessing module
- Added hard check: the observed and synthetic data format must be correct
- Added hard check: synthetic and observations files must match
- Added function to check that all required adjoint sources exist
- Critical Bugfix: Adjoint traces were incorrectly written as synthetics and not
  as calculated adjoint sources.

### #158
- Removes GitHub links in dependencies in favor of PyPi packages
- Fixes broken Pyatoa import statements
- Pyadjoint misfit functions are now listed properly in the Pyaflowa docstring- 


### #168
This is an important PR, see PR notes for more detailed description

#### MAJOR BUGS
- Parallel written residuals file could sometimes have bad formatting
- Residuals were appended to the same file during line search causing incorrect
	misfit to be calculated
- Trial model (m_try) was never exposed to line search meaning model was 
	not updated

#### Bugfixes
- Residuals files now written per event, iteration and step 
- Max step length now defined correctly
- Export residuals flag now works properly
- Synthetic traces exported at each evaluation


### #176
Support for HPC system Wisteria 

- Makes the location of the `submit` and `run` scripts a System Class variable 
	which can be overwritten. 
- Adds `system.Fujitsu` as a generalized System class for HPCs using the Fujitsu
	workload manager
- Adds `system.Wisteria` for System interactions with HPC Wisteria
- Custom run and submit scripts for Wisteria hardcoded to RSCGRP RC58, not 
	generalized at the moment

