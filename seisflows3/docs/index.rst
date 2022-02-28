.. SeisFlows documentation master file, created by
   sphinx-quickstart on Mon Oct 26 16:35:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SeisFlows3
===============

SeisFlows is a Python-based waveform inversion package used to tackle the
problems of full waveform inversion, seismic migration, and adjoint tomography.

With a growing user base in academia and industry, this package has been used
for production scale inversions, some with over a billion model parameters, for
research problems related to oil and gas exploration, earthquake seismology, and
general nonlinear optimization problems.

SeisFlows3 is a fork of SeisFlows and the current main development branch.
Changes to the package include migration to Python3, updates to source-code
readability through improved doc strings and comments, stronger adherance to
PEP-8 standards, and additional functionalities not present in the original
package.

Throughout the documentation we may use the names SeisFlows and SeisFlows3
interchangeably to refer to this package. Any reference to the original
SeisFlows will be made explicitely. Due to incompatabilities between
Python 2 and 3, backwards compatability between SeisFlows and SeisFlows3 is not
supported.


Installation
=================

Successful applications of SeisFlows3 will typically require direct editing of
source code. For this reason, SeisFlows3 should be installed directly via the
package library using pip, using the -e flag to ensure that SeisFlows3 is
installed in development mode, allowing source code changes to be immediately
acccessible to Python.

The `devel` branch houses the most up-to-date codebase. We recommend installing
SeisFlows3 within a virtual environment (e.g., Conda) to preserve your root
environment.

.. code:: bash

   $ conda create -n seisflows3 python=3.7
   $ conda activate seisflows3
   $ git clone --branch devel  https://github.com/bch0w/seisflows3/
   $ cd seisflows
   $ pip install -e .

.. note::
   At an undetermined future date we will migrate this package to the
   more permanent org page: https://github.com/seisflows. Please check back
   to this page to be alerted of this transition.

Requirements
-------------

In most production-scale workflows, SeisFlows3 must be run on a cluster, or
high performance computing system. However, serially run example problems
making use of 2D solvers like SPECFEM2D are available for small problems and
workflow tutorials.

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Introduction to SeisFlows

   overview
   background

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: How to use this package

   usage
   structure
   standards
   tips_and_tricks

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Examples and Tutorials

   tutorials

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Additional Material

   citeme
   development
   devplan
   seisflows_api

