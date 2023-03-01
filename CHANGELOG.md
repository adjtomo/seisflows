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
