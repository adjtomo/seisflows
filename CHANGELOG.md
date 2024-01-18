# SeisFlows Changelog

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


## v2.3.0
A collection of bugfixes and feature improvement (big thanks to @evcano for major PR #168)

- Improved code organization for readability 
- Allow acoustic domain and pressure seismograms in SPECFEM2D (#164)
- Improved workflow feature `export_residuals` (#168)
- Allow workflow toggling of SPECFEM parameter `SAVE_FORWARD_ARRAYS`
- System support for HPC Wisteria (U. Tokyo) (#176)
	- Makes the location of the `submit` and `run` scripts a System Class 
	  variable which can be overwritten. 
	- Adds `system.Fujitsu` as a generalized System class for HPCs using the 
	  Fujitsu workload manager
	- Adds `system.Wisteria` for System interactions with HPC Wisteria
	- Custom run and submit scripts for Wisteria hardcoded to RSCGRP RC58, not 
	  generalized at the moment

### Bugfixes
- Broken symlinks for Cluster runscript locations when SeisFlows is 
  installed via Pip (#162), added Manifest.in file to fix this.
- Max workers in concurrent jobs was set with the wrong variable (#160), replaced iwth correct variable
- Bracketing line search allows non-passing models to pass line search under certain criteria, fixed.
- Solver directory intialization was too general when searching for 
  source files (#169), leading to an incorrect number of source files being detected, process now only looks for sources files that start with source prefix (e.g., CMTSOLUTION, FORCESOLUTION)
- Residuals files were sometimes written incorrectly (#168) which lead to broken file reading. Residuals files are now written out at a more granular level to avoid multiple processes trying to write to one file
- Residuals files were appended to by all processes during line search, leading to incorrect misfit values (#168); fixed by above bugfix.
- Updated trial model was not exposed to line search after step count > 1   (#168), fixed
- Incorrect step lengths defined when resuming workflow (#168), fixed
- Exported synthetic traces were overwritten each iteration (#168). They are now tagged by iteration + step count to avoid overwriting
- output_optim.txt was not writing the correct misfit values for each 
  iteration, fixed

