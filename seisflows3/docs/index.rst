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

Throughout the documentation we use the name SeisFlows to refer to this package
unless explicitely stated otherwise. Due to incompatabilities between
Python 2 and 3, backwards compatability between SeisFlows and SeisFlows3 is not
guaranteed.


Installation
=================

Successful applications of SeisFlows will require direct editing of source code.
For this reason, SeisFlows should be installed directly via the package library
using pip.

The -e flag in the pip install command ensures that SeisFlows is
installed in development mode, meaning changes to the source code are
immediately acccessible to the Python interpreter.

.. code:: bash

   $ git clone https://github.com/seisflows/seisflows.git
   $ cd seisflows
   $ pip install -e .

Requirements
-------------

In most production scale workflows, SeisFlows must be run on a cluster, or
high performance computing system.

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
   seisflows_api

