
1. Download SeisFlows
---------------------

To install Seisflows, first clone the repository::

    git clone https://github.com/PrincetonUniversity/seisflows


Then set environment variables. If using bash, add the following lines to ``.bash_profile`` or ``.bashrc`` (or modify accordingly, if you are using a different shell)::

    export PATH=$PATH:/path/to/seisflows/scripts
    export PYTHONPATH=$PYTHONPATH:/path/to/seisflows


2. Load example
---------------

After installing SeisFlows, users with accounts on ``tiger.princeton.edu`` can type ``sfexamples`` at the command line to interactively choose an example. Once an example is chosen, a ``parameters.py`` file and ``paths.py`` file will be copied to the current directory.


3. Run example
--------------

To submit a job, type ``sfrun`` from within the directory containing ``parameters.py`` and ``paths.py`` files. If the ``serial`` system option is specified in ``parameters.py``, the job will begin executing immediately. If the ``tiger`` system configuration are specified, then the job will run under SLURM when sufficient resources become available. Once the job starts running, status information will be displayed to the terminal or to the file ``output.log``.

After trying the example once to make sure everything is working, users can explore different inversion settings by modifying entries in ``parameters.py``.  Entries in ``paths.py`` correspond to actual locations on the tigress filesystem and cannot be changed.

If you experience problems, first check that you are using the latest version of SeisFlows. To update, type `git pull` within the  directory containing SeisFlows source code.
