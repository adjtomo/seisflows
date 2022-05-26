Getting started
=================================

Assuming installation detailed on the `home page <index.html>`__ has
completed successfully, you can start playing around with SeisFlows.

To look at the command line help tool, you can type:

.. parsed-literal::

    seisflows -h

Each sub-argument has it's own help message to further explain what it does.

.. parsed-literal::

    seisflows par -h

For more information on the SeisFlows command line tool, see the
`command line tool <command_line_tool.html>`__ docs page.

Running an example problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To get further than the command line help messages, SeisFlows requires an
external numerical solver to operate.

We recommend running example problem 1, which will download and compile
`SPECFEM2D <https://geodynamics.org/cig/software/specfem2d/>`__, and then run a
simple inversion workflow.

.. parsed-literal::

    seisflows examples run 1

Interested in what's going on under the hood of this example problem? See the
`Specfem2D workstation example <specfem2d_example.html>`__ docs page.

Understanding SeisFlows
~~~~~~~~~~~~~~~~~~~~~~~~~

If the example problem has completed successfully, the directory that you're in
should resemble a standard SeisFlows working directory.

    *  Decipher the SeisFlows parameter file at the
       `parameter file <parameter_file.html>`__ docs page.
    *  Understand the working directory structure at the
       `working directory <working_directory.html>`__ docs page.
    *  For an exploration of the code structure, have a look at the
       `API <autoapi/index.html>`__ or the
       `GitHub repo <https://github.com/bch0w/seisflows>`__.
    *  Need to expand SeisFlows for your own research? Learn how to
       `extend the package <extending.html>`__, e.g., to bring in new workflows
       or teach SeisFlows how to interact with various HPC compute systems.

