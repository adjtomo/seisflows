.. SeisFlows documentation master file, created by
   sphinx-quickstart on Thu Oct 16 14:47:24 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SEISFLOWS
=========

Background
----------
The most sophisticated waveform inversion packages are proprietary codes developed by oil and gas companies and geophysical service providers.  Usually, such packages are maintained by software engineering teams and not available to independent researchers.

Outside of industry, a number of open source waveform inversion packages have been developed, mainly in geophysics but also in nondestructive testing and other areas. Many early waveform inversion packages used simple frequency-domain formulations, which were well-suited for research but not readily scalable.  

More recent packages such as SeisFlows use Python for data processing tasks in combination with parallel compiled code for wave simulation.  This approach combines the ease of use of modern scientific Python and with the efficiency and scalability of modern acoustic and elastic time-domain solvers.

With SeisFlows, wave simulations must be performed using an external software package such as SPECFEM2D or SPECFEM3D. The ability to interface with external solvers provides flexibility, and the choice of SPECFEM as a default option gives access to cutting-edge solver capabilities. However, the need for an external package imposes some additional demands on the user as described `here <http://seisflows.readthedocs.io/en/latest/instructions_remote.html#creating-your-own-examples>`_.

SeisFlows provides an automated nonlinear optimization workflow, with the option to carry out multiple model upates without stopping or to stop between updates for quality control checks.


Examples: Available For Download
--------------------------------

We have prepared a 2D waveform inversion example that is inexpensive enough to run on almost any laptop, desktop, or cluster. To run this simple checkboard inversion, see these step-by-step :doc:`instructions <instructions_remote>`. 

Some additional examples are available for download.  Please review the instructions for the 2D checkerboard test case to get a sense for how to run these other inversions.

Some 2D examples based on the Marmousi model are available `here <http://tigress-web.princeton.edu/~rmodrak/2dElastic>`_.

A 3D Cartesian checkerboard example is available `here <http://tigress-web.princeton.edu/~rmodrak/3dElastic>`_.

A 3D global 1-chunk example is available `here <http://tigress-web.princeton.edu/~rmodrak/ExamplesGlobal>`_. Please note, the compressed archive for this example is very large (> 0.5 GB). [No longer available because of file size.]

At a minimimum, one processer is required for the 2D Marmousi examples, 16 processors are required for the 3D Cartesian example, and 64 processors are required for the global 1-chunk example.  See `here <http://seisflows.readthedocs.org/en/latest/usage/usage.html#system-configuration>`_ for more information about running inversions in parallel.

Note: File hosting services are provided by my alma mater.  The download server may become temporarily unavailable due to system maintenance or permanently unavailable due to expiration of my account.


Examples: Available Locally
---------------------------

Users with accounts on "tiger.princeton.edu" can run the following inversions without having to download files or recompile executables.

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

   usage/contents.rst

