SeisFlows
=========

SeisFlows is an open source seismic inversion package that

- Delivers a customizable waveform inversion workflow

- Provides a framework for research in regional, global, and exploration seismology

Examples and usage information are available at [readthedocs.org](http://seisflows.readthedocs.org/en/latest/).  Probably the easiest way to learn more is to follow these checkerboard example [instructions](http://seisflows.readthedocs.io/en/latest/instructions_remote.html).

The most sophisticated waveform inversion packages are proprietary codes developed by oil and gas companies. Usually, such software is not available to independent researchers.  Outside of industry, a number of open source waveform inversion packages have been developed.  SeisFlows differs from previous open-source software in providing a high-level inversion workflow emphasizing both flexibility and HPC portability.

With SeisFlows, the inversion task is abstracted into six components: `solver`, `system`, `nonlinear optimization`, `data preprocessing`, `postprocessing`, and `workflow`.  This design is informed by hands-on experience with different HPC environments and research applications. The source code is structured in a modular way based on these six categories, and users are offered various choices in each one.  To see the choices available for each category, simply browse the source code.  The inversion itself is executed by `seisflows/workflow/inversion.py`, which may be good place to start.

Wave simulations must be performed using an external software package. The ability to interface with external solvers provides flexibility, and the choice of SPECFEM2D/3D as default options gives optional GPU acceleration and other useful capabilities. However, the need for an external package imposes some additional demands on the user.  Setting up your own inversions using SPECFEM2D/3D can be time consuming because it requires familiarity with SPECFEM2D/3D's idiosyncratic meshing procedure and binary file formats. Alternatively, some users have interfaced with others external solvers, but this can also be time consuming.


References
----------
If you find this package useful, please cite:

`Ryan Modrak, Jeroen Tromp; Seismic waveform inversion best practices: regional, global and exploration test cases, Geophysical Journal International, Volume 206, Issue 3, 1 September 2016, Pages 1864–1889, https://doi.org/10.1093/gji/ggw202`

Another manuscript is currently under review in Computers and Geosciences.


See also
--------
The following extension packages are not currently documented, but may still give a sense for the type of research possible within this framework:

- https://github.com/rmodrak/seisflows-research
- https://github.com/rmodrak/seisflows-multiparameter
- https://github.com/rmodrak/seisflows-hpc



[![Build Status](https://travis-ci.org/rmodrak/seisflows.svg?branch=master)](https://travis-ci.org/rmodrak/seisflows)


