Extending the Codebase
======================

.. note::
    Page Under Construction
  
SeisFlows works on the object oriented programming principal of inheritance.
See the `inheritance <inheritance.html>`__ page for background. This allows
Users to extend the package to work with other systems and solvers, or modify
its behavior to suit the problem at hand.

Overriding SeisFlows with your own subclass
-------------------------------------------

If the existing modules of SeisFlows do not suit your needs, you may
need to override them with your own sub classes. A common example of
when you might need to do this is to override the system subclasses to
tailor SeisFlows to your specific cluster.

In this example SeisFlows already contains a Slurm system module, and
has overriding sub-classes for Slurm-derived systems including Chinook
(University of Alaska Fairbanks), Maui (New Zealand eScience
Infrastucture) and Frontera (Texas Advanced Computing Center). Any new
Slurm systems will need write new subclasses to take advantage of the
existing structure.

Below we provide some examples of editing existing SeisFlows modules for
modified performance.


Writing your own Base class from scratch
----------------------------------------

Less likely, but also still a possibility, is that you will have to
write your own Base class which defines foundational capabilities. One
example of this is extending SeisFlows to interact with other numerical
solvers. Currently SeisFlows is set up to work with
SPECFEM2D/3D/3D_GLOBE. To allow Seisflows to work with other solvers,
one would need to write a new Base class that defines SeisFlows’
interaction behavior with this solver.


Setting up 'Cluster' System sub-classes
-----------------------------------------

In this section, we outline how various Cluster-derived System sub-classes
are created. Although the tasks performed here will not translate one-to-one
for other systems, we hope they serve as a roadmap for setting up SeisFlows on
other HPC systems.

Because many HPC systems have their own unique caveats, it is standard practice
to trial and error our way into a working sub-class.

.. note::
    It is assumed that SeisFlows has been installed with the
    ``pip install -e .`` command (develop mode). This ensures that any changes
    to the codebase are immediately propagated into the SeisFlows executable.

University of Alaska Fairbanks 'Chinook'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Chinook is University of Alaska Fairbank's high performance computer. At the
time of writing (Aug. 24, 2022), it is running CentOS 6.10 and the SLURM
workload manager.

Because Chinook is a SLURM-based system, we can derive our 'Chinook' class
based on the SLURM system class. In seisflows/system we create a new
Python file called 'chinook.py', which has a base template of:

.. code:: python3

    from seisflows.system.slurm import Slurm

    class Chinook(Slurm):
        """
        Chinook System - UAF SLURM-based HPC

        Parameters
        ==========

        Paths
        =====

        ***
        """
        __doc__ == Slurm.__doc__ + __doc__

        def __init__(self):
            """Chinook intitiation"""
            super().__init__()


The above code snippet is boilerplate code for any new SeisFlows class. We can
see that `Chinook` inherits base functionalities from `Slurm`, it defines some
boilerplate class docstring and appends it's docstring to the `Slurm` class
(required by ``seisflows configure``), and inherits the SLURM startup
procedure defined in __init__.

.. note::
    Because of an outdated operating system, we cannot run the
    SeisFlows directly from Chinook's Conda installation. We therefore have to
    containerize the codebase and run SeisFlows through 'Singularity/Apptainer'.

We can now set up a TestFlow workflow to see how the Chinook implementation
differs from the bog-standard SLURM impementation.

.. code:: bash

    mkdir ${CENTER1}/scratch  # CENTER1 is the Chinook work filesystem
    cd ${CENTER1}/scratch
    seisflows setup  # create the parameters.yaml file
    # Set the core SeisFlows modules
    seisflows par workflow test_flow
    seisflows par system chinook
    seisflows par solver null
    seisflows par preprocess null
    seisflows par optimize null

    seisflows configure  # set up the remainder of the parameter file

