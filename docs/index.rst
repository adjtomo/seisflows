.. SeisFlows documentation master file, created by
   sphinx-quickstart on Thu Oct 16 14:47:24 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

What is SeisFlows?
==================

SeisFlows is an open source seismic inversion package that:

- delivers a complete, customizable adjoint tomography and full waveform inversion workflow

- provides a framework for research in regional, global, and exploration tomography


Installation
============

To install Seisflows, first clone the repository::

    git clone https://github.com/PrincetonUniversity/seisflows


Then set environment variables. If using bash, add the following lines to ``.bash_profile`` or ``.bashrc`` (or modify accordingly, if you are using a different shell)::

    export PATH=$PATH:/path/to/seisflows/scripts
    export PYTHONPATH=$PYTHONPATH:/path/to/seisflows


Getting Started
===============

Perhaps the easiest way to get started is to run a few examples.  Users with accounts on ``tiger.princeton.edu`` can run a number of ready-to-go test cases following these :doc:`instructions <tiger>`.  

Other users can download and run an example following the alternate :doc:`instructions <instructions>`.


Documentation
=============

.. toctree::
   :maxdepth: 1

   manual/contents.rst
   ref/modules.rst

