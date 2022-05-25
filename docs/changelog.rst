Change Log
===============

The following list documents changes made from the original ``SeisFlows``
codebase (Last updated April 18, 2022).

Major
------
* 'seisflows' command line tool which now acts as the main entry point into
  package. Replaces the old sf* scripts and adds functionality to inspect
  the package and an active workflow.
* Overhauled the parameter file, merging the old parameters.py and paths.py into
  a single parameters.yaml file.\
* Overhauled path and parameter input system to be stricter and more informative
* Created a seisflows-super directory directly inside the package, replacing
  the previous style of maintaining separate repos for custom modules.
* Package-wide implementation of Python Logger module, replacing old-style
  print statements. Improved verbosity and inheritance tracking.
* Renamed source code namespace from "seisflows" to "seisflows"
* Updated all syntax to Python3.7
* Modernized setup.py file which sets the console script command line tool
  during install and mandates package requirements


Moderate
--------
* Renamed and restructured SeisFlows working directory structure.
    - output.logs -> logs
    - output.stats -> stats
    - output.optim -> stats/line_search.txt
* ``stop_after`` and ``resume_from`` parameters allow a user to (as suggested)
  STOP a workflow after a certain workflow function (e.g., forward simulations)
  and RESUME a workflow from a certain function.
* Added template parameter file, base class and super class files to make it
  easier to further develop the package
* Removed seisflows.plugins.optimize and shifted ALL functionality into
  Optimize module to centralize optimization functionality.
* Mandated file extensions on output models and text files to make it easier
  to distinguish file types. For example, models USED to be saved as e.g.,
  m_new (Numpy .npy files). They are now saved as e.g., m_new.npy
* Removed the "_sm" and "_lg" system distinctions.
* Dynamic parameter file generation based on the new SeisFlowsPathsParameters
  class. That is, based on the user-chosen modules, a custom parameter file
  which is ONLY relevant to the chosen modules is written out with docstrings,
  type hinting, and default values.
* Overhauled parameter checking check() functions in ALL modules.
  Check functions no longer allowed to set or overwrite parameters, rather
  they inform users when parameters are set incorrectly.
* Removed any "default" module classes in liue of using the "base" classes.
* Renamed seisflows.tools.tools to seisflows.tools.wrappers to avoid
  folder and file name repetition
* Renamed seisflows.tools.seismic to seisflows.tools.specfem to make it
  clearer that these tools are for SPECFEM-related operations.
* Overhauled the seisflows.tools.specfem.getpar() and setpar() functions for
  getting and setting values in the SPECFEM Par_file
* Separated system mid-tier classes into "Cluster" and "Workstation" to
  differentiate working on HPC systems


Minor
------
* reST format docstrings for all classes and functions.
* Comments added throughout codebase to explain naming, logic etc.
* Replaced all string formatters with f-strings if possible, .format() otherwise
* Replaced literal path concatenations ('path' + '/' + 'to') with os.path.join()
* Removed all unncessary function abstractions (e.g., function loadnpy() simply
  called np.load(); replaced with with np.load())
* Stronger adherance to PEP-8 including: CamelCase classes, snake_case modules
  and functions.
* Overhauled the seisflows.config.custom_import function to adhere imports to
  PEP-8
* Updated subprocess calling to use currently accepted API run() function, as
  opposed to previous implementation using check_call() and check_ouput() etc.
* Merged the system.run() and system.run_single() into a single run() function,
  avoids what was essentially copy-pasted code with small tweaks.
* Replaced individual Writer classes which were attributes of optimization and
  preprocess modules. These provided an unncessary layer of abstraction from
  simple file writing.
* Enforced package-wide constants at the top of seisflows.config
* Added __init__() functions to most of the modules to define any
  instance-dependent variables, which previously were not explained.


