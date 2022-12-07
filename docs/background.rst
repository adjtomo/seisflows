Background
==========
Learn how the SeisFlows source code is structured and the design choices that
went into building the package.

--------------------------------

Modularity
-------------
SeisFlows is built on five core modules which each take care of a separate
functionality. Most modules are independent, only relying on their own set of
paths and parameters.

The exception to this is the `Workflow` module, which calls and combines the
other modules within its own functions to solve the task at hand.

Within each module, a **base** Class defines the general functionality.
Additional 'child' classes inherit and modify the 'base' class.
The SeisFlows Modules are:

* **Workflow**: Controls the order and execution of tasks performed during a
  workflow, e.g., first create working directory, *then* run simulations.
  (Base: *seisflows.workflow.forward.Forward*)
* **System**: Compute system interface used to run SeisFlows on different
  computer architectures. A consistent internal structure makes it relatively
  seamless to switch between workstation and HPC workload manager
  implementations.
  (Base: *seisflows.system.workstation.Workstation*)
* **Solver**: Interface and wrapper for external numerical solvers used to
  generate models, synthetic waveforms, and gradients.
  (Base: *seisflows.solver.specfem.Specfem*)
* **Preprocessing**: Signal processing operations performed on time series,
  as well as adjoint source generation and misfit windowing.
  (Base: *seisflows.preprocess.default.Default*)
* **Optimization**: Nonlinear optimization algorithms used to find
  the minimum of a given objective function.
  (Base: *seisflows.optimize.gradient.Gradient*)


--------------------------------


.. include:: inheritance.rst


--------------------------------


An Example of Inheritance in SeisFlows
----------------------------------------

Here we look at how we define the System module ``Chinook``, which allows
SeisFlows to interface with the University of Alaska supercomputer, Chinook.

- **System.Workstation**: ``Workstation`` is the **Base** system class. It
  defines required paths and parameters as well as general functions for
  running SeisFlows on a workstation/laptop.

- **System.Cluster**: The ``Cluster`` class inherits properties from
  ``Workstation`` and additionally modifies some of the call behavior so that
  SeisFlows can run on HPCs. For example, ``Cluster`` defines new parameters:
  **WALLTIME** (dictates the length of time a workflow is submitted for) and
  **TASKTIME** (defines the expected time required for any single simulation).

- **System.Slurm** The ``Slurm`` system inherits from ``Cluster``, defining
  even **more** specific functionality, such as checking of the job queue
  using the ``sacct`` command.

.. note::

    Note that those using other workload managers (e.g., PBS, LSF), will need to
    branch off here and define a new class that inherits from ``Cluster``.

- **System.Chinook**: Individual Slurm systems differ in how Users interact
  with the workload manager. ``Chinook`` defines the calls structure for
  working on the Slurm-based University of Alaska Fairbanks HPC, named Chinook.

  ``Chinook`` inherits **all** of the attributes of ``Workstation``, ``Cluster``
  and ``Slurm``, and adds it's own specifics such as the available partitions
  on the cluster.

- **SeisFlows System**: SeisFlows abstracts all of this behavior behind a
  generalized ``System`` module. Therefore, calling something like
  ``system.submit()`` will bring call on ``Chinook.submit()``, which may
  inherit from any of the parent classes.

  All SeisFlows modules are similarly structured, with a **Base** classes defining
  standard behavior, and Child classes extending and modifying such behavior.

