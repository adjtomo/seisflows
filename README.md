SeisFlows
=========

SeisFlows is an open source seismic inversion package that

- Delivers a complete, customizable waveform inversion workflow

- Provides a framework for research in regional, global, and exploration seismology

Examples and usage guidelines are available online at [readthedocs.org](http://seisflows.readthedocs.org/en/latest/).

SeisFlows has an idiosyncratic workflow-based design informed by familiarity with different research problems (oil and gas exploration, earthqake tomography, medical imaging) and experience in different computing environments, including a variety of high-performance computing clusters.

With SeisFlows, wave simulations must be performed using an external software package such as SPECFEM2D or SPECFEM3D.  The ability to interface with external solvers ensures flexibility, and the choice of SPECFEM as a default option gives access to cutting-edge meshing and hardware accelaration capabilities.  At the same time, the need for an external solver creates additional work for the user.  

To carry out an inversion, one must become familiar not only with SeisFlows, but also with a separate solver package.  If you decide to try SeisFlows, it should be relatively straightforward to run the checkerboard and Marmousi examples available at [readthedocs.org](http://seisflows.readthedocs.org/en/latest/).  To move beyond these examples, one would need to become familiar with how to set up simulations with SPECFEM, and in paricular, with how to create models in SPECFEM's idionsyncratic binary format.


References
----------
If you use this package in your research, please cite:

`Ryan Modrak, Jeroen Tromp; Seismic waveform inversion best practices: regional, global and exploration test cases, Geophysical Journal International, Volume 206, Issue 3, 1 September 2016, Pages 1864–1889, https://doi.org/10.1093/gji/ggw202`

Another manuscript is currently under review in Computers and Geosciences.


Relation to other packages
--------------------------
SeisFlows is one of a few dozen waveform inversion packages.  Most of these are sophisticated, proprietary codes developed by oil and gas companies and geophysical service providers.  Generally, such packages are maintained by professional software engineering teams and not available to independent researchers.

Outside of industry, a number of open source packages have been developed, mainly in geophysics but also in nondestructive testing and other areas.  The terms _waveform inverison_, _full waveform inversion_, _adjoint tomography_, and _finite frequeny tomography_ are all largely synonymous, and associated software shares many similarities.

Many early waveform inversion packages used simple frequency-domain formulations, which were well-suited for research but not readily scalable beyond inexpensive 2D problems.  More recent packages such as SeisFlows use Python for data processing tasks in combination with parallel compiled code for wave simulation.  This approach combines the ease of use of modern scientific Python and with the efficiency and scalability of modern time-domain acoustic and elastic wave-equation solvers.

Unlike some packages that require extensive human intervention, SeisFlows provides an automated nonlinear optimization workflow, with the option to carry out multiple model upates without stopping or to stop between updates for quality control checks.


See also
--------
Many specialized inversion strategies have been implemented by overloading SeisFlows classes.  The following extension packages give a sense for the type of research possible within the framework:

- https://github.com/rmodrak/seisflows-research
- https://github.com/rmodrak/seisflows-multiparameter


The following extension package illustrates how SeisFlows can be tailored to specific high-performance computing environments and how  fault tolerance can be implemented:

- https://github.com/rmodrak/seisflows-hpc



[![Build Status](https://travis-ci.org/rmodrak/seisflows.svg?branch=master)](https://travis-ci.org/rmodrak/seisflows)
