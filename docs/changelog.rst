Change Log
===============

The following list documents changes made since v2.0.0 
(Last updated Aug. 18, 2022).

Major
------
* Reworks the entire parameter system to take advantage of Python's native init
  and internal attributes, rather than setting global parameters in sys.modules.
* Removes reliance on globally accessible sys.modules, which was used so that all
  modules had access to one another, but also prevented standalone imports of modules
  and efficient unit testing. This has been replaced with explicit I/O for modules to
  talk to each other, and shared parameters that are instantiated by a single parameter
  file.
* Sys.modules was also used to checkpoint workflows, as all class instances were saved
  as pickle files that could be re-loaded. However this required a specific startup
  procedure anytime a workflow had to be instantiated. This has been replaced by an
  explicit checkpointing system which takes advantage of simpler text files which save
  only relevant checkpoint information.
* Implemented a unit/integration test structure using Pytest, for all modules and utilities.
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
* Support added for SLURM-based TACC Frontera system
* All modules now have their own setup() and check() functions which help
  establish a workflow before anything needs to be submitted to the system.
* System run functions no longer require global pickling into sys modules, but
  instead simply pickle the required functions and their keyword arguments. This
  greatly simplifies the run capabilities on external systems, and removes the need
  to define specific pickling functions (copyreg) within the package.
* Solver has been separated from all model manipulations. Replaced by a new
  standalone Model class which completely characterizes a Specfem model and can
  read/write models as numpy vectors or as Specfem-readable binary files.
* Postprocess module has been scrapped and its functionalities absorbed by the
  Solver and Workflow modules. The reasoning behind this was that the postprocess
  really had no agency of its own, and was just a wrapper for solver functionality.
  It may be reintroduced in the future if we find need for very complex postprocessing
  steps.
* Optimization test suite now runs Rosenbrock problem (from legacy code) for all
  optimization schema and tests line search capabilities.
* Optimization also checkpoints itself and the line search during a workflow, allowing
  Users to restart workflows from any point effortlessly and with minimal loss of compute time.
  Pyatoa preprocessing renamed to 'Pyaflowa'. All of the old Pyaflowa class functions
  (previously contained in the Pyatoa package) are now exposed in the SeisFlows Pyaflowa
  preprocessing class. This is so that users do not necessarily have to look at the Pyatoa
  source code to understand what's going on in SeisFlows
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
* Added a TestFlow workflow which 'live' tests System abilities directly on
  login/compute nodes of a cluster
* Reworked some of the internal system.Cluster architecture to take advantage
  of inheritance - Systems now specify headers for their run and submit calls.
  e.g., slurm-based systems specify their own SBATCH commands but all use the 
  same call structure after the SBATCH header


Minor
------
* Removed a lot of agency from individual modules (i.e., writing their own files,
  manipulating other modules) and gave it mainly to the Workflow module. This is so
  that when users need to look at where changes are disk are happening, they only have
  to look at the Workflow module, rather than chasing through various other
  modules/sub-modules.
* Rearranged directory structure to be more easily navigable, with less top-directory scripts.
* SeisFlows now relies more heavily on inheritance and super() functions throughout,
  whereas previously it was used quite sparingly.
* Workflow module now relies more heavily on inheritance. Forward workflows build out most
  of the functionality of the forward solver, and the Migration and Inversion workflows build
  on top of this.
* Reworked docstring structure to aid in 'seisflows configure' command which now builds new
  parameter files directly from init and class docstrings
* Gutted most of the old internal start up procedures such as loading parameters and modules.
  The new start up procedure is now more Pythonic and efficient.
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
* New Docs page to describe how to transition from a 2D workstation example 
  to a 2D cluster example.

