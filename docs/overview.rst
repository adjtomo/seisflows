Overview
===================
``SeisFlows`` is a Python-based package for automating waveform inversion.
With a user base in both academia and industry, ``SeisFlows`` has been used
for:
- production scale inversions (some with over a billion model parameters),
- research problems related to oil and gas exploration,
- regional-scale earthquake-based adjoint tomography problems,
- general nonlinear optimization problems.

Why use SeisFlows?
~~~~~~~~~~~~~~~~~~

SeisFlows is meant to reduce the *"human-time"* cost of performing a seismic
inversion. It replaces the typical hodgepodge of Bash, Python or Matlab
scripts Users inevitably write to run a seismic inversion on a supercomputer.

In addition, it is open-source, community-driven and has been used in published
studies.

How does it work?
~~~~~~~~~~~~~~~~~~
With only a parameter file, a User can type a single command on an HPC terminal
and launch an automated inversion workflow onto their cluster. Under the hood,
a main job acts on behalf of the User: submitting simulations, monitoring the
job queue, managing the filesystem, and generating results.

``SeisFlows`` is built upon the object-oriented programming concept of
`inheritance <background.html>`__. Architecturally, SeisFlows hides
case-by-base details behind a generalized namespace such that a workflow can
remain consistent whether applied to a 2D acoustic problem on a Linux
workstation, or to a continental-scale adjoint tomography problem run on a
high performance computer.

Thanks to inheritance, if desired functionality is missing in the  package,
new users can overload default classes without having to rewrite entire
portions of the code base, allowing SeisFlows to maintain a rigid structure
while still being flexible in its application to various scales and problems.


What does SeisFlows not do?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``SeisFlows`` is **not** a numerical solver, **not** a waveform processing tool,
and **not** a visualization tool.
Knowledge of interfacing with external numerical solvers, as well as waveform
preprocessing tools, is paramount to successful application of SeisFlows.


How do I use SeisFlows?
~~~~~~~~~~~~~~~~~~~~~~~

Check out the `Getting Started <getting_started.html>`__ page to learn how to
run an example problem, understand a SeisFlows working directory and set up
SeisFlows for your own research problems.

-------------------------

What is SeisFlows Legacy?
--------------------------

SeisFlows Legacy is the name given to the initial Python2 release of the
codebase. The current version of SeisFlows has seen a number of
backwards-incompatible changes from this initial version.

For Users who used SeisFlows prior to these updates, note that major
backwards-incompatible changes from the legacy codebase include:

-  complete shift to Python>=3.7 source code, abandoning Python2 support
-  reworking of internal code architecture to be more 'Pythonic' and modular
-  richer source code emphasizing readability and PEP-8 standards
-  a new command line tool for improved package control
-  redesigned, dynamically-generated parameter file

See the `change log <changelog.html>`__ for point-by-point changes from the
original codebase. The `Legacy SeisFlows codebase
<https://github.com/adjtomo/seisflows/releases/tag/v1.1.0>`__ can still be
accessed along with its `documentation <https://github.com/adjtomo/seisflows/
tree/46a9604d8907bc5fbf2ddc11e16914a870e6db77/docs>`__, but is no longer
supported by the developers.