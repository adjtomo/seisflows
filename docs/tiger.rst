
Loading examples
----------------

After installing SeisFlows, users with accounts on ``tiger.princeton.edu`` can type ``sfsetup`` at the command line to interactively choose an example. Once an example is chosen, a ``parameters.py`` file and ``paths.py`` file are created in the current directory.


Running examples
----------------

To submit a job, type ``sfrun`` from within the directory containing ``parameters.py`` and ``paths.py`` files. If the 'serial' system option is specified in ``parameters.py``, the job will begin executing immediately. If ``pbs`` or ``slurm`` system configurations are specified, then the job will run when resources become available. Once the job starts running, status information will be displayed either to the terminal or to the file ``output.log``.

After trying the example once to make sure everything is working, users can explore different inversion settings by modifying the values in ``parameters.py``. Values in ``paths.py``, on the other hand, correspond to locations on the tigress filesystem and cannot be changed.


