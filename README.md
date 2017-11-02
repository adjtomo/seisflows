SeisFlows
=========

SeisFlows is an open source seismic inversion package that

- Delivers a customizable waveform inversion workflow

- Provides a framework for research in regional, global, and exploration seismology

Examples and usage information are available at [readthedocs.org](http://seisflows.readthedocs.org/en/latest/).  Probably the easiest way to learn more is to follow these checkerboard example [instructions](http://seisflows.readthedocs.io/en/latest/instructions_remote.html).

References
----------
If find this package useful, please cite:

`Ryan Modrak, Jeroen Tromp; Seismic waveform inversion best practices: regional, global and exploration test cases, Geophysical Journal International, Volume 206, Issue 3, 1 September 2016, Pages 1864–1889, https://doi.org/10.1093/gji/ggw202`

Another manuscript is currently under review in Computers and Geosciences.


Design
------
The inversion task is abstracted into six components: solver, system, nonlinear optimization, data preprocessing, image postprocessing, and workflow.  This design is informed by hands-on experience with many different high-performance computing environments and research applications (including oil and gas exploration, earthquake tomography, and ultrasound imaging).  

To see the choices available for each component, simply browse the source code and note the modules available in each of the six directories.  The inversion itself is executed by `inversion.main` in `seisflows/workflow/inversion.py`, which may be another good place to browse.



See also
--------
The following extension packages are not currently documented, but may still give a sense for the type of research possible within this framework:

- https://github.com/rmodrak/seisflows-research
- https://github.com/rmodrak/seisflows-multiparameter
- https://github.com/rmodrak/seisflows-hpc



[![Build Status](https://travis-ci.org/rmodrak/seisflows.svg?branch=master)](https://travis-ci.org/rmodrak/seisflows)
