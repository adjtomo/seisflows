.. SeisFlows documentation master file, created by
   sphinx-quickstart on Thu Oct 16 14:47:24 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

What is SeisFlows?
==================

SeisFlows is an open source seismic inversion package with three main goals:

- deliver a complete, customizable adjoint tomography and full waveform inversion workflow

- provide a framework for doing research in regional, global, and exploration tomography

- avoid code duplication across 2D, 3D, and global solver packages


Installation
============

To install Seisflows, first clone the repository::

    git clone github.com/PrincetonUniversity/seisflows


Then set environment variables. If using bash, add the following lines to ``.bash_profile`` or ``.bashrc``::

    export PATH=$PATH:/path/to/seisflows/scripts
    export PYTHONPATH=$PYTHONPATH:/path/to/seisflows


Getting Started
===============

Perhaps the easiest way to get started is to run a few examples.  Users with accounts on ``tiger.princeton.edu`` can start directly with these :doc:`instructions <tiger>`.  Other users must await the arrival of an ftp download server. 


Documentation
=============

.. toctree::
   :maxdepth: 1

   manual/contents.rst
   ref/modules.rst

