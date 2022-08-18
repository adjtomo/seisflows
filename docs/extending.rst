Extending SeisFlows
===================

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

