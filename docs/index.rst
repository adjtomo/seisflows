.. SeisFlows documentation master file, created by
   sphinx-quickstart on Thu Oct 16 14:47:24 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

What is SeisFlows?
==================

SeisFlows is an open source seismic inversion package that

- delivers a complete, customizable waveform inversion workflow

- provides a framework for research in regional, global, and exploration seismology

Perhaps the easiest way to learn more is to run a one of the following examples.  


Examples: Available For Download
--------------------------------

We have prepared a 2D waveform inversion example that is inexpensive enough to run comfortably on almost any laptop, desktop, or cluster. To run this simple checkboard inversion, see these step-by-step :doc:`instructions <instructions_remote>`. 

Some additional examples are available for download.  Please review the instructions for the 2D checkerboard test case to get a sense for how to run these other inversions. Aside from differences in computational cost and under-the-hood `solver <http://seisflows.readthedocs.org/en/latest/manual/manual.html#solver-configuration>`_ and `system <http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration>`_ details, each of these inversion examples runs in basically the same way.

Some 2D examples based on the Marmousi model are available `here <http://tigress-web.princeton.edu/~rmodrak/2dElastic>`_.

A 3D Cartesian checkerboard example is available `here <http://tigress-web.princeton.edu/~rmodrak/3dElastic>`_.

A 3D global 1-chunk example is available `here <http://tigress-web.princeton.edu/~rmodrak/Examples3dGlobe>`_. Please note, the compressed archive for the 3D global example (1 GB) is much larger than for all the other examples.

At a minimimum, one processer is required for the 2D Marmousi examples, 16 processors are required for the 3D Cartesian example, and 64 processors are required for the global 1-chunk example.  See `here <http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration>`_ for more information about running inversions in parallel.



Examples: Available Locally
---------------------------

Users with accounts on "tiger.princeton.edu" can run the following inversions without having to download files are recompile executables.

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


Documentation
=============

.. toctree::
   :maxdepth: 1

   manual/contents.rst

