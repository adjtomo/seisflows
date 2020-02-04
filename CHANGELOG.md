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
