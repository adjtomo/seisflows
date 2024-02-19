.. image:: images/sf_globe_banner_alpha.png
    :align: center

------------------

SeisFlows --- Automated Inversion Software
==========================================

`SeisFlows <https://github.com/adjtomo/seisflows>`__  is an open-source,
Python-based waveform inversion package meant to automate the problems of full
waveform inversion, seismic migration, and adjoint tomography on high
performance computers.

SeisFlows is hosted on `GitHub <https://github.com/adjtomo/seisflows>`__ as
part of the `adjTomo organization <https://github.com/adjtomo/>`__.

---------------------------------

Quickstart
~~~~~~~~~~~

- Check out the `Overview page <overview.html>`__ to learn about SeisFlows and
  how to use it.
- Found a bug or want to request a new feature? Feel free to
  `open a GitHub Issue! <https://github.com/adjtomo/seisflows/issues>`__
- Want to talk about SeisFlows?
  `Check in on the discussions page. <https://github.com/orgs/adjtomo/discussions>`__
- Contributions are encouranged and welcome! Have a look at the 
  `adjTomo Contributor's Guide <https://pyatoa.readthedocs.io/en/latest/contributing.html>`__  to see how you can contribute (hosted on Pyatoa's docs).


---------------------------------

Installation
~~~~~~~~~~~~

To install SeisFlows and its dependencies, we recommend installing within
a Conda environment to not affect your root environment. The 
`master <https://github.com/adjtomo/seisflows/>`__ branch contains the most
stable release and is the preferred installation branch. The 
`devel <https://github.com/adjtomo/seisflows/tree/devel>`__ branch houses 
the most up-to-date codebase, but is likely to be **unstable**.

.. note::

    SeisFlows is installed using the Pip ``-e`` which enables development
    mode, where source code changes are immediately acccessible to Python (i.e.,
    you do **not** need to re-install SeisFlows when updating source code).

.. code:: bash

   git clone https://github.com/adjtomo/seisflows.git
   cd seisflows
   conda env create -f environment.yml
   conda activate seisflows

Alternative: update an existing environment
```````````````````````````````````````````

The above installation instructions create a **new** Conda environment named 
`seisflows`. If you have an *existing* Conda environment and want to install 
SeisFlows and it's dependencies there, you can run the following:

.. code:: bash

    git clone https://github.com/adjtomo/seisflows.git             
    cd seisflows 
    conda activate <your environment>
    conda env update --name <your environment> --file environment.yml

---------------------------------


Cite SeisFlows
~~~~~~~~~~~~~~~~~~

If you use SeisFlows for your work, please cite the following papers:
`Chow et al. (2020) <https://doi.org/10.1093/gji/ggaa381>`__ and
`Modrak et al. (2018)
<https://www.sciencedirect.com/science/article/pii/S0098300417300316>`__

    Chow, B., Kaneko, Y., Tape, C., Modrak, R., & Townend, J. (2020).
    *An automated workflow for adjoint tomography --- waveform misfits and
    synthetic inversions for the North Island, New Zealand.*
    Geophysical Journal International, 223(3), 1461-1480.

    Modrak, R. T., Borisov, D., Lefebvre, M., & Tromp, J. (2018).
    *SeisFlows — Flexible waveform inversion software.*
    Computers & geosciences, 115, 88-95.

Publications which have used (and cited) SeisFlows can be
found on `Google Scholar <https://scholar.google.com/scholar?cites=9435477750683593672&as_sdt=405&sciodt=0,2&hl=en>`__.


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Introduction

   overview
   getting_started
   background

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
   running_on_chinook

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: How To's

   tips_and_tricks
   cluster_setup
   containers

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Development

   extending
   changelog
   code_dev_plan

