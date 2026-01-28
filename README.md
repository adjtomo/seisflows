SeisFlows 
==========

[![Documentation Status](https://readthedocs.org/projects/seisflows/badge/?version=devel)](https://seisflows.readthedocs.io/en/devel/?badge=devel)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

SeisFlows is an open-source, Python-based waveform inversion package that tackles the problems of seismic wavefield simulation, full waveform inversion, seismic migration, and adjoint tomography on a variety of compute systems, from laptops to supercomputers. Seisflows is built on top of external numerical solvers to make it easier for users to run large, complicated workflows with lots of moving parts.

The **design principles** of SeisFlows include: 
1. Driven through: command line Interface + dynamically generated parameter file.
2. Inheritance-based design for generalized workflows and plug-and-play functionality.
3. Automated workflows with checkpointing, pre-written job scheduler interactions, and verbose logging.

Interested?
------------

- Please see [Read the Docs](https://seisflows.readthedocs.io) for install instructions, examples, and API documentation.

- SeisFlows is housed package under the [adjTomo organization](https://github.com/adjtomo) which contains tools for computational seismology. SeisFlows has strong ties to the [SPECFEM](https://specfem.org/) community.


- If you find any issues, have questions, or would like to join the community, please feel free to open up a [GitHub Issue](https://github.com/adjtomo/seisflows/issues) or [start a discussion](https://github.com/orgs/adjtomo/discussions). 


References
----------
If you use this package in your own research, please cite the following papers:

- Bryant Chow, Yoshihiro Kaneko, Carl Tape, Ryan Modrak, John Townend, *An automated workflow for adjoint tomography     —waveform misfits and synthetic inversions for the North Island, New Zealand*, Geophysical Journal International, Volume 223, Issue 3, December 2020, Pages 1461–1480, https://doi.org/10.1093/gji/ggaa381

- Ryan Modrak, Dmitry Borisov, Matthieu Lefebvre, Jeroen Tromp; *SeisFlows—Flexible waveform inversion software*, Computers & Geosciences, Volume 115, June 2018, Pages 88-95, https://doi.org/10.1016/j.cageo.2018.02.004

