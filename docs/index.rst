.. SeisFlows documentation master file, created by
   sphinx-quickstart on Thu Oct 16 14:47:24 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

What is SeisFlows?
==================

SeisFlows is an open source seismic inversion package that

- delivers a complete, customizable waveform inversion workflow

- provides a framework for research in regional, global, and exploration seismology


Getting Started
===============

Perhaps the easiest way to get started is to run a one of the following examples.  


Available Locally
-----------------

Users with accounts on tiger.princeton.edu can run the following ready-to-go inversions.

*2D Regional and Global*

- North America

- Southern California

- Global

- Deep Earth

*2D Near Surface*

- Marmousi offshore

- Marmousi onshore

- overthrust offshore

- overthrust onshore

- BP anticline

- BP salt diapir

*3D Cartesian*

- checkerboard

*3D Global*

- mideast

See these :doc:`instructions <instructions_local>` for running inversions on our local cluster.


Available For Download
----------------------

Users without accounts on tiger.princeton.edu can download and run the following inversions.

*2D*

- checkerboard

See these :doc:`instructions <instructions_remote>` for downloading and running inversions.  

The examples above can be run on most any laptop, desktop, or cluster. For more information about running on various systems, see `here <http://seisflows.readthedocs.org/en/latest/instructions_remote.html#run-checkerboard-test-in-parallel>`_ and `here <http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration>`_.

We are working on making more examples available for download.  Once you are familiar, however, with the working of both the package and the solver (for the above examples, SPECFEM2D), it can be straightforward to set up inversions of your own.



Documentation
=============

.. toctree::
   :maxdepth: 1

   manual/contents.rst

