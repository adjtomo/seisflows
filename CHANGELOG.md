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

