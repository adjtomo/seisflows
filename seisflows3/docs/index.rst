.. SeisFlows documentation master file, created by
   sphinx-quickstart on Mon Oct 26 16:35:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


SeisFlows3
===========

`SeisFlows3 <https://github.com/bch0w/seisflows3>`__  is a hard fork of 
`SeisFlows <https://github.com/rmodrak/seisflows>`__, a Python-based
waveform inversion package used to tackle the problems of full waveform
inversion, seismic migration, and adjoint tomography. SeisFlows3 constitutes
the current main development branch of SeisFlows.

With a user base in both academia and industry, ``SeisFlows`` has been used
for: production scale inversions (some with over a billion model parameters),
research problems related to oil and gas exploration, regional-scale
earthquake-based adjoint tomography problems, and general nonlinear
optimization problems.


.. note::
   SeisFlows\ **3** is written in Python\ **3**. Major backwards-incompatible
   changes from the original SeisFlows codebase include:

      * complete shift to Python3.7 source code, abandoning Python2 support
      * richer source code emphasizing readability and standards
      * a new command line tool for improved package control
      * redesigned, dynamically-generated parameter file
      * native integration with the waveform misfit quantification tool: `Pyatoa <https://github.com/bch0w/pyatoa>`__

   See the `change log <changelog.html>`__ for point-by-point changes from the
   original codebase.


Throughout the documentation we use the names
`seisflows`, `SeisFlows`, `SeisFlows3` etc. interchangeably to refer to
this package (SeisFlows3). Any reference to the original SeisFlows will be
noted explicitely.

---------------------------------

Installation
=================

.. note::
   At a yet-to-be determined date we will migrate this package to the
   more permanent org page: https://github.com/seisflows.  Please check back
   to this page to be alerted of this transition.

Successful applications of SeisFlows3 will typically require direct editing of
source code. For this reason, SeisFlows3 should be installed directly via the
package library using Pip. The ``-e`` flag ensures that SeisFlows3 is
installed in development mode, allowing source code changes to be immediately
acccessible to Python.

The `devel` branch houses the most up-to-date codebase. We recommend installing
SeisFlows3 within a virtual environment (e.g., Conda) to preserve your root
environment.

.. code:: bash

   $ conda create -n seisflows3 python=3.7
   $ conda activate seisflows3
   $ git clone --branch devel https://github.com/bch0w/seisflows3.git
   $ cd seisflows
   $ pip install -e .

--------------------------------

Requirements
=============

In most production-scale workflows, SeisFlows3 must be run on a cluster, or
high performance computing system. However, serially run example problems
making use of 2D solvers like SPECFEM2D are available for small problems and
workflow tutorials.

SeisFlows3 + Pyatoa
--------------------
To include the waveform measurement capabilities of Pyatoa, you must install
separately. See the `Pyatoa documentation
<https://pyatoa.readthedocs.io/en/latest/>`__ for the most up to date
install instructions.

.. code:: bash

   $ cd ..
   $ git clone --branch devel https://github.com/bch0w/pyatoa.git
   $ cd pyatoa
   $ pip install .


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Introduction

   overview
   background
   citeme

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: How To

   usage
   structure
   standards
   tips_and_tricks

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Examples

   specfem2d_example
    
.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Development

   changelog
   development
   code_dev_plan
#   seisflows_api

