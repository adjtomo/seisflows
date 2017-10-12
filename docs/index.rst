.. SeisFlows documentation master file, created by
   sphinx-quickstart on Thu Oct 16 14:47:24 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SeisFlows is an open source seismic inversion package that

    Delivers a complete, customizable waveform inversion workflow

    Provides a framework for research in regional, global, and exploration seismology

SeisFlows has an unique, somewhat idiosyncratic workflow-based design informed by familiarity with different research problems (including oil and gas exploration, earthqake tomography, medical imaging) and experience working in different computing environments (including various high-performance computing clusters).

Probably the easiest way to learn more is to run one of the examples below.  


Relation to other packages
--------------------------
The most sophisticated waveform inversion packages are expensive proprietary codes developed by oil and gas companies and geophysical service providers.  Generally, such packages are maintained by professional software engineering teams and not available to independent researchers.

Outside of industry, a number of open source packages have been developed, mainly in geophysics but also in nondestructive testing and other areas. Many early waveform inversion packages used simple frequency-domain formulations, which were well-suited for research but not readily scalable beyond 2D problems.  More recent packages such as SeisFlows use Python for data processing tasks in combination with parallel compiled code for wave simulation.  This approach combines the ease of use of modern scientific Python and with the efficiency and scalability of modern acoustic and elastic wave-equation solvers.

SeisFlows provides an automated nonlinear optimization workflow, with the option to carry out multiple model upates without stopping or to stop between updates for quality control checks.

With SeisFlows, wave simulations must be performed using an external software package such as SPECFEM2D or SPECFEM3D. The ability to interface with external solvers provides flexibility, and the choice of SPECFEM as a default option gives access to cutting-edge solver capabilities. However, the need for an external package imposes some additional demands on the user as described `here <http://seisflows.readthedocs.io/en/latest/instructions_remote.html#creating-your-own-examples>`_.


Examples: Available For Download
--------------------------------

We have prepared a 2D waveform inversion example that is inexpensive enough to run on almost any laptop, desktop, or cluster. To run this simple checkboard inversion, see these step-by-step :doc:`instructions <instructions_remote>`. 

Some additional examples are available for download.  Please review the instructions for the 2D checkerboard test case to get a sense for how to run these other inversions.

Some 2D examples based on the Marmousi model are available `here <http://tigress-web.princeton.edu/~rmodrak/2dElastic>`_.

A 3D Cartesian checkerboard example is available `here <http://tigress-web.princeton.edu/~rmodrak/3dElastic>`_.

A 3D global 1-chunk example is available `here <http://tigress-web.princeton.edu/~rmodrak/ExamplesGlobal>`_. Please note, the compressed archive for this example is very large (> 0.5 GB).

At a minimimum, one processer is required for the 2D Marmousi examples, 16 processors are required for the 3D Cartesian example, and 64 processors are required for the global 1-chunk example.  See `here <http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration>`_ for more information about running inversions in parallel.



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

   manual/contents.rst

