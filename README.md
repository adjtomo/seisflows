SeisFlows
=========

SeisFlows is an open source seismic inversion package that

- Delivers a complete, customizable waveform inversion workflow

- Provides a framework for research in regional, global, and exploration seismology

Examples and usage guidelines are available online at [readthedocs.org](http://seisflows.readthedocs.org/en/latest/).

SeisFlows has an idiosyncratic workflow-based design informed by familiarity with diverse research problems (oil and gas exploration, earthqake tomography, medical imaging) and experience running in diverse Unix environments, including a wide variety of high-performance computing clusters.

With SeisFlows, wave simulations must be performed using an external software package such as SPECFEM2D or SPECFEM3D.  The ability to interface with external solvers ensures flexibility.  The choice of SPECFEM2D/3D for the default solver interfaces creates an opportunity to take advantage of cutting-edge meshing and hardware accelaration capabilities.

At the same time, the need for an external solver creates additional work for the user.  To carry out an inversion, one must become familiar not only with SeisFlows, but also with a separate solver package.  If you decide to try SeisFlows, it should be fairly straightforward to run the checkerboard and Marmousi examples available at [readthedocs.org](http://seisflows.readthedocs.org/en/latest/).  To move beyond these examples, one would need to become familiar with how to set up simulations with SPECFEM2D/3D, and in paricular, with how to create models in SPECFEM2D/3D's idionsyncratic binary format.


References
----------
If you use this package in your research, please cite:

`Ryan Modrak, Jeroen Tromp; Seismic waveform inversion best practices: regional, global and exploration test cases, Geophysical Journal International, Volume 206, Issue 3, 1 September 2016, Pages 1864–1889, https://doi.org/10.1093/gji/ggw202`

Another manuscript is currently under review in Computers and Geosciences.


Relation to other packages
--------------------------
SeisFlows is one of perhaps a few dozen waveform inversion packages.  The most sophisticated of these are expensive, proprietary codes developed by leading oil and gas companies and geophysical service providers.  Generally, such packages are maintained by professional software engineering teams and not available to independent researchers.

Outside of industry, several open source packages have been developed, mainly in geophysics but also in nondestructive testing and other areas.  The terms _waveform inverison_, _full waveform inversion_, _adjoint tomography_, and _finite frequeny tomography_ are all more or less synonymous, and associated software packages share many similarities.

Many early waveform inversion packages used simple frequency-domain formulations, which were suitable for research but not scalable beyond inexpensive 2D problems.  More recent packages such as SeisFlows use Python for inexpensive processing tasks in combination with parallel compiled code for wave simulation.  This approach combines the ease of use of modern scientific Python and with the efficiency and scalability of modern time-domain acoustic and elastic wave-equation solvers.  Unlike some packages, SeisFlows provides an fully automated nonlinear optimization workflow, which is useful for prototyping inversion methods.


See also
--------
Many specialized inversion strategies have been implemented by overloading SeisFlows classes.  The following extension packages give a sense for the type of research possible within the framework:

- https://github.com/rmodrak/seisflows-research
- https://github.com/rmodrak/seisflows-multiparameter


The following extension package illustrates how SeisFlows can be tailored to specific high-performance computing environments and how  fault tolerance can be implemented:

- https://github.com/rmodrak/seisflows-hpc



[![Build Status](https://travis-ci.org/rmodrak/seisflows.svg?branch=master)](https://travis-ci.org/rmodrak/seisflows)
