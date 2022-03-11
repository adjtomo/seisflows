Background
==========================
Understanding what SeisFlows is meant for and how it works requires a level of
background understanding. Here we provide brief explanations aimed to providing
a general overview of topics required to understand SeisFlows. Where applicable
we provide links and references to further, more in-depth material.


Inheritance in Object Oriented Programming
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Inheritance is one of the main concepts of object oriented programming. It
allows users to create new classes that share some, or all, of the
characteristics of an existing class. This allows new users to build upon
existing work, reducing the amount of efforts required to tailor existing
classes to their specific needs.

* Baseclass: Class that defines its own attributes and does not inherit
* Superclass: A class being inherited from (Parent)
* Subclass: A class that inherits some or all of its attributes (Child)

A class can be both a super and a subclass if there are multiple levels of
inheritance.

To do:
 * Create a small Jupyter notebook that shows off how inheritance and the
   super() function work in Python.
 * Links to further reading related to inheritance in Python and OOP


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

