.. SeisFlows documentation master file, created by
   sphinx-quickstart on Mon Oct 26 16:35:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: images/sf_globe_banner_alpha.png
    :align: center

------------------

`SeisFlows <https://github.com/adjtomo/seisflows>`__  is a Python-based
waveform inversion package used to tackle the problems of full waveform
inversion, seismic migration, and adjoint tomography. SeisFlows constitutes
the current main development branch of SeisFlows.

With a user base in both academia and industry, ``SeisFlows`` has been used
for: production scale inversions (some with over a billion model parameters),
research problems related to oil and gas exploration, regional-scale
earthquake-based adjoint tomography problems, and general nonlinear
optimization problems.


.. note::
   Major backwards-incompatible changes from the legacy codebase include:

      -  > complete shift to Python3.7 source code, abandoning Python2 support
      -  > reworked internal architecture to be more 'Pythonic'
      -  > richer source code emphasizing readability and standards
      -  > a new command line tool for improved package control
      -  > redesigned, dynamically-generated parameter file
      -  > native integration with the waveform misfit quantification tool:
         `Pyatoa <https://github.com/adjtomo/pyatoa>`__

   See the `change log <changelog.html>`__ for point-by-point changes from the
   original codebase.


.. warning::
    This docs page is currently under active development and may have missing
    information or unfinished pages. We are doing our best to complete it in a
    timely manner.

---------------------------------

Installation
=================

Successful applications of SeisFlows will typically require direct editing of
source code. For this reason, SeisFlows should be installed directly via the
package library using Pip. The ``-e`` flag ensures that SeisFlows is
installed in development mode, allowing source code changes to be immediately
acccessible to Python.

The `devel` branch houses the most up-to-date codebase. We recommend installing
SeisFlows within a virtual environment (e.g., Conda) to not affect your root
environment.

.. code:: bash

   $ conda create -n seisflows python=3.10
   $ conda activate seisflows
   $ git clone --branch devel https://github.com/adjtomo/seisflows.git
   $ cd seisflows
   $ conda install --file requirements.txt
   $ pip install -e .

--------------------------------

Requirements
=============

In most production-scale workflows, SeisFlows must be run on a cluster, or
high performance computing system. However, serially run example problems
making use of 2D solvers like SPECFEM2D are available.

SeisFlows + Pyatoa
--------------------
To include the waveform measurement capabilities of Pyatoa, you must install
separately. See the `Pyatoa documentation
<https://pyatoa.readthedocs.io/en/latest/#installation>`__ for the most up to date
install instructions.


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Introduction

   overview
   start_here
   citeme

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Structure

   command_line_tool
   parameter_file
   working_directory

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Examples

   specfem2d_example
   example_on_cluster

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: How To's

   containers
   extending

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Development

   background
   design
   changelog
   code_dev_plan

