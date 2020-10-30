## 1.1.0
* Updated all syntax to Python3, changed import statements to reflect some 
  package namespace changes from Python2
* Docstrings and parameter type and descriptions for all classes and functions
* Explanation docstrings for each module, and explanations on which modules are
  base modules and which modules inheret
* Overhauled the parameters.yaml file to include type and descriptions for each
  parameter, as well as section descriptions. Also includes custom parameters
  for user set Parameters
* Completely removed the parameters.py and paths.py files in lieu of a 
  parameters.yaml file, which now contains both parameters and paths. 
  These can just be set to null or commented out and the submit script will 
  purge them from the internally used parameters file.
* Removed the bash based sf* scripts (e.g. sfsubmit) in lieu of a python
  based seisflows script that provides a bit more control over the functionality
  - Replaces functionality of sfsubmit, sfresume, sfclean 
  - Includes restart functionality and a status check 
  - Included bash calls (e.g. sfsubmit) to call python script for backward 
    compatability
  - All future interaction with seisflows can go through this script
* Reformatted code to adhere loosely to PEP-8
  - All classes are now CamelCase, modules, functions are now snake_case
  - Config custom_import function can now search for snake_case modules and will
    look for the corresponding class as CamelCase. Future custom classes will
    need to adhere to this naming schema
* Added solver and workflows that work with Pyatoa without the need for
  subprocess calls
* Pedantic but replaced all calls to 'os' with full package calls, e.g. 
  os.path.join() rather than from os.path import join; join(). This makes for 
  easier identification of paths being called within the code
* Replaced literal path concatenations ('path' + '/' + 'to') with os.path.join()
* Replaced some print statement warnings with actual warnings package calls

## 1.2.0
* Completely overhauled the parameter setting functionality by mandating 
  parameters be set using a SeisFlowsPathsParameters class, which requires
  docstrings, type and default values to be set. 
* Overhauled parameter settings now allows for dynamic generation of the 
  parameter file, meaning the user doesnt have to hunt down individual 
  parameters anymore, and with the mandate of docstrings, users can 
  automatically determine type and use of each parameter that is set
* Changed name of seisflows.config.config() to init_seisflows() for a more 
  descriptive name of the function and to avoid naming confusion
* Modern setup.py file that provides command line entry points to the package
* Command line seisflows tool that allows easier interaction with the package
  and provides things like init and config functions that streamline the 
  initiation of a seisflows workflow
* seisflows-super directory now present in the main package repository, 
  which removes the need to 'export PYTHONPATH' to another directory and keeps
  all development in a single main repo
* Totally new documentation page that will include tutorials, examples and API
