Package Structure
====================
SeisFlows generates its own internal directory structure and saves the state
of the workflow at each critical step. This docs page </structure>
explains how one navigates the SeisFlows directories to efficiently run your
inversion.

The SeisFlows class
~~~~~~~~~~~~~~~~~~~
All command-line commands are passed through the SeisFlows class, which exists
in seisflows3/scripts/seisflows.py. If you ever ask yourself, "I wonder what
commands are being executed when I enter: seisflows submit, you need only go to
the SeisFlows class and look at the 'submit' function. This will tell you
exactly what commands are being executed.

SeisFlows Modules
~~~~~~~~~~~~~~~~~~~~~

Working Directory Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scratch Directory
~~~~~~~~~~~~~~~~~~~

Outputs
~~~~~~~~~~~

Pickled Workflow State
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
