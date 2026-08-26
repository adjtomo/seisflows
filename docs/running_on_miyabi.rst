Running on JCAHPC Miyabi
=========================

This page introduces the SeisFlows ``Miyabi`` system interface, and how to
configure a SeisFlows parameter file to run on
`Miyabi <https://miyabi.jcahpc.jp/>`__, the supercomputer operated by the
Joint Center for Advanced High Performance Computing (JCAHPC), a
collaboration between the University of Tsukuba and the University of
Tokyo.

.. note::

    This page was written against the Miyabi User's Guide v1.8 (Fujitsu
    Limited, 2025). Contact JCAHPC support or consult the latest User's
    Guide on the `User Portal <https://miyabi-www.jcahpc.jp/>`__ if
    anything here appears to be out of date.

0. Overview
~~~~~~~~~~~~

Miyabi runs the **PBS Professional** job management system, so SeisFlows
interacts with it through the ``Miyabi`` system class (``Workstation`` →
``Cluster`` → ``Pbs`` → ``Miyabi``), which sits alongside the existing
SLURM-based (``Chinook``) and Fujitsu/PJM-based (``Wisteria``) interfaces.

Miyabi is split into two independently-scheduled subsystems that do **not**
share a login node or a CPU architecture. Which one you use depends on
whether your SPECFEM build targets GPU or CPU execution:

.. list-table::
   :header-rows: 1
   :widths: 15 40 40

   * -
     - Miyabi-G
     - Miyabi-C
   * - Node hardware
     - 1x NVIDIA GH200 Grace-Hopper Superchip (Arm Neoverse V2 CPU, 72
       cores + 1x H100 GPU)
     - 2x Intel Xeon MAX 9480 CPU (56 cores each, 112 cores/node)
   * - Architecture
     - aarch64 (Login-G)
     - x86_64 (Login-C)
   * - Node count
     - 1,120
     - 190

You must log in to, and build your solver on, the login node matching the
subsystem you intend to run on (``miyabi-g.jcahpc.jp`` or
``miyabi-c.jcahpc.jp``) -- a binary compiled on Login-G will not run on
Miyabi-C's compute nodes, and vice versa.

.. warning::

    The ``Miyabi`` system module has not yet been run against a live
    allocation. Please run `TestFlow <cluster_setup.html#testflow>`__ first
    (Section 3 below) to validate job submission/monitoring on your
    account before attempting a full simulation or inversion, and please
    `open a GitHub Issue <https://github.com/adjtomo/seisflows/issues>`__
    with any quirks you run into so the interface can be improved.

1. Prerequisites
~~~~~~~~~~~~~~~~~

- SSH access to Miyabi (public key + OTP, set up through the JCAHPC
  `User Portal <https://miyabi-www.jcahpc.jp/>`__)
- A project code with available "tokens" (Miyabi's core-hour equivalent).
  Check your allocation on the login node with:

  .. code:: bash

      show_token

- SeisFlows and SPECFEM2D/3D/3D_GLOBE installed and compiled *on the
  matching login node* (Login-G for Miyabi-G, Login-C for Miyabi-C). See
  the `main installation instructions <index.html#installation>`__. Miyabi
  loads Conda as a module by default (check ``module avail``), so
  environment creation should work as normal once the module is loaded.

.. note::

    Miyabi's job management system does not support multibyte (non-ASCII)
    characters anywhere in a job submission -- job names, project codes,
    directory names, etc. The ``Miyabi`` system module enforces this for
    its own parameters (``title``, ``group``, ``queue``) in ``check()``,
    but this restriction also applies more broadly (e.g., to your working
    directory path).

2. Set the Conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like Wisteria, Miyabi's compute nodes do **not** inherit the login node's
Conda environment, so a plain ``qsub``-submitted job will not have access
to the Python environment SeisFlows is installed in. Set the ``conda_env``
parameter to your Conda environment name (or full path), and SeisFlows will
route job submission/execution through the shared
``runscripts/conda_activate-miyabi`` wrapper script, which activates that
environment on the compute node before running SeisFlows:

.. code:: bash

    seisflows par conda_env seisflows
    # or, e.g., a full path:
    seisflows par conda_env /work/<group>/<user>/conda/envs/seisflows

If ``conda_env`` is left unset, jobs are submitted directly with no
wrapping -- only appropriate if Conda is made available to compute node
jobs some other way (e.g., inherited via ``-V``, if your site allows it).

.. note::

    This mechanism (and this module's defaults more generally) is
    currently only set up/tested for Miyabi-G (GPU). Running on Miyabi-C
    should work by selecting a ``-c`` ``queue``, but has not been
    exercised.

3. Validate with TestFlow
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before running real simulations, confirm SeisFlows can submit array jobs,
monitor the PBS queue, and correctly detect job failures on your account
using `TestFlow <cluster_setup.html#testflow>`__:

.. code:: bash

    mkdir testflow_miyabi && cd testflow_miyabi
    seisflows init

    seisflows par workflow test_flow
    seisflows par system miyabi
    seisflows par solver null
    seisflows par preprocess null
    seisflows par optimize null

    seisflows configure

4. Configure the ``Miyabi`` parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running ``seisflows configure`` will populate your parameter file with all
``Miyabi``-specific parameters (in addition to the ones inherited from
``Pbs``/``Cluster``/``Workstation``). The two you must set are:

.. code:: bash

    # Project code used to charge tokens (PBS '-W group_list=')
    seisflows par group <your project code>

    # Node-occupied use queue, e.g. one of:
    #   debug-g, short-g, regular-g (Miyabi-G, aarch64 + GH200)
    #   debug-c, short-c, regular-c (Miyabi-C, x86_64)
    seisflows par queue debug-g

Other parameters worth double-checking against your job size (see the
``Miyabi`` and ``Pbs`` class docstrings, or ``seisflows par -h queue`` /
``seisflows par -h walltime`` etc. for the full description of each):

.. code:: bash

    seisflows par ntask 2       # number of events/tasks
    seisflows par nproc 72      # cores per task, matches SPECFEM's nproc
    seisflows par tasktime 30   # walltime (minutes) per spawned job
    seisflows par walltime 60   # walltime (minutes) for the main job
    seisflows par mpiexec mpirun  # default; matches Miyabi-G's NVIDIA HPC SDK

.. note::

    Only "node-occupied use" queues are currently supported (whole compute
    nodes). Miyabi-G's fractional-GPU "MIG use" queues
    (``debug-mig``/``short-mig``/``regular-mig``) are not implemented, since
    SPECFEM-style multi-node MPI workflows are expected to use whole nodes.

5. Submit
~~~~~~~~~~

As with any other SeisFlows system, the master job knows how to:

- submit jobs (using ``qsub``)
- monitor the queue (using ``qstat``, falling back to ``qstat -H`` and
  ``tracejob`` to determine whether a finished job succeeded or failed,
  since PBS's own job STATUS does not distinguish the two)
- stop the workflow gracefully if a job fails and cannot be retried

.. code:: bash

    seisflows submit

Monitor progress the same way as on any other system: the main log is
written to ``sflog.txt``, and each spawned job writes its own log file to
``logs/``.

.. warning::

    Do **not** pass ``--direct`` to ``seisflows submit``/``restart`` on
    Miyabi. That flag is only meaningful for systems (like Wisteria) whose
    ``submit()`` supports an in-process fallback mode; ``Pbs``/``Miyabi``
    do not implement it, and the CLI will print a clean "does not accept
    argument `direct`" error rather than submitting anything.

6. Troubleshooting
~~~~~~~~~~~~~~~~~~~~

- **"PBS 'queue' must match one of [...]"** -- the ``queue`` parameter must
  be one of Miyabi's node-occupied use queues (see Section 4 above, or
  the ``Miyabi`` class docstring for the current list and their node/
  walltime limits).
- **"qsub: invalid option -- '-'" (or similar, mentioning single
  characters)** -- this means something is appending extra command line
  arguments after the submitted script name; PBS's ``qsub`` does not
  reliably forward those to the script (see the note in ``system.Pbs``'s
  module docstring). If you see this, you're likely running a modified/
  older version of the ``Pbs``/``Miyabi`` classes -- update to the current
  version, which passes all such information via qsub's ``-v`` option
  instead.
- **"CondaError: Run 'conda init' before 'conda activate'"** -- the batch
  job's shell never sourced Conda's shell hook. This is handled
  automatically by ``runscripts/conda_activate-miyabi`` (via
  ``eval "$(conda shell.bash hook)"``); if you still see this, confirm
  ``conda_env`` is set (Section 2 above) so that wrapper is actually being
  used, and that plain ``conda`` (not just ``python``) is on ``PATH`` once
  its module is loaded.
- **Jobs rejected at submission** -- check ``show_token`` for remaining
  allocation, and confirm your ``group`` project code is correct.
- **Multibyte character errors** -- see the note in Section 1; rename your
  working directory/job title to plain ASCII.
- For anything else, check the individual job log files in ``logs/``, and
  compare against ``qstat -H -t <job_id>`` / ``tracejob <job_id>`` run
  directly on the login node.

Have a look at `Extending SeisFlows <extending.html>`__ if you need to
further customize the ``Miyabi`` interface (e.g., to add MIG-instance
queue support, or a custom ``mpiexec`` per subsystem).
