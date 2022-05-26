Background
==========================
Understanding what SeisFlows is meant for and how it works requires a level of
background understanding. Here we provide brief explanations aimed to providing
a general overview of topics required to understand SeisFlows. Where applicable
we provide links and references to further, more in-depth material.

.. include:: inheritance.rst

Inheritance in SeisFlows
----------------------------

SeisFlows employs inheritance almost exactly like the example above. Take the
``System`` module for example. The Base class defines required paths (such as
the location of the OUTPUT directory on disk. The definition of this path is
system-independent (it's just a path). The setup() function is similarly
system-independent, as it just creates the directory structure in the working
directory.

One Child of our ``System`` Baseclass is the ``Cluster`` module, which defines
attributes for a module that interacts with high performance compute systems.
The ``Cluster`` module defines more specific parameters, such as a WALLTIME,
dictating the length of time a workflow is submitted for.

For those working on a ``Slurm`` system, we create a child of ``Cluster``.
``Slurm`` defines even **more** specific functionality, such as job status
checking using the ``sacct`` command.

And finally, we can have a specific Child which caters to the unique cluster
we're working on. For example, ``Chinook`` defines the calls structure for
working on the University of Alaska Fairbanks HPC, named Chinook. ``Chinook``
inherits **all** of the attributes of ``Base``, ``Cluster`` and ``Slurm``,
reducing the amount of code repitition, keeping required structure consistent
for all of the Parent classes, while still providing the maximum amount of
flexibility of working on a specific compute system.

SeisFlows abstracts all of this behavior behind a generalized ``System``
module. Therefore, calling something ``system.submit()`` will bring in all
of the inherited properites of the above mentioned system superclasses.

All SeisFlows modules are similarly structured, with ``Base`` classes defining
standard behavior, and Child classes extending and modifying these behavior.


..
    !!! TO DO
    Numerical Waveform Simulations using SPECFEM
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SeisFlows acts as a Python-based interface for SPECFEM2D, SPECFEM3D Cartesian
    and SPECFEM3D Globe. However, to access the functionalities contained in these
    solvers, the user must separately retrieve the SPECFEM source code and build
    the associated binary files. SPECFEM has an extensive documentation page which
    should guide your through these tasks. Once these are generated, SeisFlows only
    needs to know the location of the SPECFEM DATA directory, containing the
    associated parameter, station and source files, as well as the SPECFEM BIN
    directory, which contains the executable binary files.

    To do:
     * Small code snippet for downloading/compiling SPECFEM from GitHub
     * Include relevant links to SPECFEM documentation page


    Theory of Full Waveform Inversion
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Full waveform inversion and adjoint tomography are seismic imaging techniques
    used to improve characterizations of subsurface Earth structure. Many
    publications are avaialble explaining the the theory and application of seismic
    inversion. These tomographic techniques are an amalgamation of various
    scientific fields, including signal processing, numerical approximation,
    nonlinear optimization, and inverse theory.

    To do:
     * Include Fig. 5 from Chow et al. (workflow figure) and use this to broadly
       explain how adjoint tomography works
     * Explain how a normal inversion works, explain how thrifty inversion differs.
     * Relevant links to papers by Fichter, Tape, Komatitsch & Tromp, Modrak to
       explain the required theories behind FWI.

