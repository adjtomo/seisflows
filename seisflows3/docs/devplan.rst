Development Plan
===========================
A short- and long-term to-do list for SeisFlows3 development. This is a community document meant to highlight what
the SeisFlows3 developers deem as priorities for the package, allowing community members to target their contributions
more effectively.

The following plan sets the roadmap for the one year mark (11/22). The main objective for this dev plan is a
stable release of SeisFlows3, featuring a comprehensive documentation page, extensive unit and integration testing, and
usable working examples for new users.

Expected difficulty (i.e., easy, medium, hard) of each problem is provided in bracketed terms, e.g., [easy].
This difficulty is subjective and based on required Python/coding experience, package familiarity, and time required.
Tasks are also categorized, e.g., (general), when possible.

Version 1.0 (Updated Nov. 18, 2021)
-----------------------------------
* (bookkeeping) [medium]: Migrate a stable release (v1.0) of SeisFlows 3 to github.com/seisflows

* (stable) [hard]: Implement package-wide unit testing using pytest
* (stable) [hard]: Implement integration testing using 2D example problems
* (stable) [medium]: Outline and complete documentation.
* (stable) [medium]: Complete package install using Conda.

* (stable) [hard]: Complete and test Preprocess sub-module for:

    * Default
* (stable) [hard]: Complete and test Solver sub-modules for:

    * SPECFEM2D
    * SPECFEM3D_GLOBE
* (stable) [hard]: Complete and test System sub-modules for:

    * Chinook (UAF)
    * Local
    * Oakforest-PACS (JCAHPC)

* (example) [hard]: SPECFEM2D examples following legacy SeisFlows example problems, on Local and HPC
* (example) [hard]: SPECFEM3D examples following legacy SeisFlows example problem, HPC
* (example) [easy]: Convert NZNorth example problem README to a .rst doc file.

* (updates): Parameter File and parameter system

    * [easy]: Add a T0 parameter to the solver.base class to define when the SOLVER thinks time 0 is.
    * [easy-medium] Increased verbosity throughout the entire package, make it sing.
    * [hard]: Auto-fill available functions from Workflow module in the RESUME_FROM and STOP_AFTER docstrings
    * [easy]: PAR.SLURMARGS and PAR.ENVIRONS need to be set as empty strings, cannot be NoneType
    * [easy]: Attenuation needs to be set-able in parameter file, currently hard-coded on
    * [hard]: Add a SAVE_PREPROCESS option to save all files within the scratch/preprocess dir.
    * [medium] Add an 'advanced' flag to all parameters and a flag in the seisflows.configure() command to
    * [easy]: seisflows.setup() should stamp git version in the parameter file

* (feature) [medium]: Copy SPECFEM output logs somewhere after each call to the solver, e.g., output_generate_databases.txt. Currently these are overwritten each time.
* (blue-sky) [uber-hard]: Based on the modules, generate a flowchart showing each class and function that is called throughout the workflow, in order.

