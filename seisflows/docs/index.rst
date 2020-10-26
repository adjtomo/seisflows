.. SeisFlows documentation master file, created by
   sphinx-quickstart on Mon Oct 26 16:35:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SeisFlows3
===============

SeisFlows is a Python-based waveform inversion package with a growing user base
in academia and industry. The package has been used for production scale
inversions, some with over a billion model parameters, for research problems
related to oil and gas exploration, earthquake seismology, and general
nonlinear optimization problems.

SeisFlows3 is the next iteration of SeisFlows, refactored in Python3, with
streamlined functionalities and additional capabilities that were not present
in the original SeisFlows package. Due to the incompatability of Python2 and
Python3, backwards compatability between SeisFlows3 and SeisFlows is not
guaranteed.


Installation
=================

In most cases, successful applications of SeisFlows require direct editing of
the source code. For this reason, SeisFlows should be installed directly
via the package library using pip. The -e flag in the pip install command
ensures that SeisFlows is installed in development mode.

.. code:: bash

   $ git clone https://github.com/seisflows/seisflows.git
   $ cd seisflows
   $ pip install -e .

Requirements
-------------

Explain the hardware and sfotware requirements for running SeisFlows.


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

