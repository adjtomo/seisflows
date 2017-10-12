SeisFlows
=========

SeisFlows is an open source seismic inversion package that

- Delivers a complete, customizable waveform inversion workflow

- Provides a framework for research in regional, global, and exploration seismology

Examples and usage guidelines are available online at [readthedocs.org](http://seisflows.readthedocs.org/en/latest/).

SeisFlows has an idiosyncratic workflow-based design informed by familiarity with different research problems (oil and gas exploration, earthqake tomography, medical imaging) and experience working in different computing environments (including various high-performance computing clusters).

With SeisFlows, wave simulations must be performed using an external software package such as SPECFEM2D or SPECFEM3D.  The ability to interface with external solvers ensures flexibility, and the choice of SPECFEM as a default option gives access to cutting-edge meshing and hardware accelaration capabilities.  However, the use of external package also creates additional work for the user, because to carry out an inversion, one must become familiar not only with SeisFlows, but also with a separate solver package.  


References
----------
If you use this package in your research, please cite:

`Ryan Modrak, Jeroen Tromp; Seismic waveform inversion best practices: regional, global and exploration test cases, Geophysical Journal International, Volume 206, Issue 3, 1 September 2016, Pages 1864–1889, https://doi.org/10.1093/gji/ggw202`

Another manuscript is currently under review in Computers and Geosciences.


See also
--------
Many specialized inversion strategies have been implemented by overloading SeisFlows classes.  The following extension packages give a sense for the type of research possible within the framework:

- https://github.com/rmodrak/seisflows-research
- https://github.com/rmodrak/seisflows-multiparameter


The following extension package illustrates how SeisFlows can be tailored to specific high-performance computing environments and how  fault tolerance can be implemented:

- https://github.com/rmodrak/seisflows-hpc



[![Build Status](https://travis-ci.org/rmodrak/seisflows.svg?branch=master)](https://travis-ci.org/rmodrak/seisflows)
