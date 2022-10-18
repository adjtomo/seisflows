.. image:: images/sf_globe_banner_alpha.png
    :align: center

------------------

SeisFlows Documentation
=======================

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

      -  > complete shift to Python>=3.7 source code, abandoning Python2 support
      -  > reworked internal architecture to be more 'Pythonic'
      -  > richer source code emphasizing readability and standards
      -  > a new command line tool for improved package control
      -  > redesigned, dynamically-generated parameter file
      -  > native integration with the waveform misfit quantification tool:
         `Pyatoa <https://github.com/adjtomo/pyatoa>`__

   See the `change log <changelog.html>`__ for point-by-point changes from the
   original codebase. The `Legacy SeisFlows codebase
   <https://github.com/adjtomo/seisflows/releases/tag/v1.1.0>`__ can still be 
   accessed along with its `documentation <https://github.com/adjtomo/seisflows/
   tree/46a9604d8907bc5fbf2ddc11e16914a870e6db77/docs>`__, but is no longer
   supported by the developers.
 


.. warning::
    This docs page is currently under active development and may have missing
    information or unfinished pages. We are doing our best to complete it in a
    timely manner.

---------------------------------

Installation
=================

``SeisFlows`` is installed using a combination of Pip and 
`Conda <https://conda.io/projects/conda/en/latest/index.html>`__. In the
future we aim to unify this into a single install command.

*Preamble: many packages will require the conda-forge channel.*

.. code:: bash
    
   conda config --add channels conda-forge

To install ``SeisFlows`` and it's dependencies, we recommend installing within 
a Conda environment to not affect your root environment. The `devel` branch 
houses the most up-to-date codebase.  

.. code:: bash

   conda create -n seisflows python=3.10
   conda activate seisflows
   git clone --branch devel https://github.com/adjtomo/seisflows.git
   cd seisflows
   conda install --file requirements.txt
   pip install -e .


SeisFlows requires the waveform measurement capabilities of 
`Pyatoa <https://github.com/adjtomo/pyatoa>`__, which currently must be 
installed manually to the same Conda environment.

.. code:: bash

   cd ..
   git clone --branch devel https://github.com/adjtomo/pyatoa.git
   cd pyatoa
   conda install --file requirements.txt
   pip install -e .


.. note::
    Successful applications of SeisFlows will typically require editing the
    source code. For this reason, SeisFlows is installed using the Pip ``-e`` 
    which enables development mode, where source code changes are immediately
    acccessible to Python.


.. note::
    In most production-scale workflows, SeisFlows must be run on high 
    performance computing systems running external numerical solvers like 
    SPECFEM3D. The User will need to install and compile these separately. 
    However, example problems in SeisFlows make use of, and can install, 
    SPECFEM2D.



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
   2D_example_walkthrough

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: How To's

   example_on_cluster
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

