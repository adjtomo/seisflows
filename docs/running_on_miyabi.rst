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
    anything here appears to be out of date. Code and documentation written
    by Claude Code and reviewed and modified by adjTomo dev team.

1. Overview
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


2. Prerequisites
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

3. Set the Conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like Wisteria, Miyabi's compute nodes do **not** inherit the login node's
Conda environment, so a plain ``qsub``-submitted job will not have access
to the Python environment SeisFlows is installed in. Set the ``conda_env``
parameter to your Conda environment name (or full path), and SeisFlows will
route job submission/execution through the shared
``runscripts/conda_activate-miyabi`` wrapper script, which activates that
environment on the compute node before running SeisFlows:

.. code:: bash

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

4. Validate with TestFlow
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

5. Configure the ``Miyabi`` parameters
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

.. note::

    Only "node-occupied use" queues are currently supported (whole compute
    nodes). Miyabi-G's fractional-GPU "MIG use" queues
    (``debug-mig``/``short-mig``/``regular-mig``) are not implemented, since
    SPECFEM-style multi-node MPI workflows are expected to use whole nodes.

6. Submit
~~~~~~~~~~

As with any other SeisFlows system, the master job knows how to:

- submit jobs (using ``qsub``)
- monitor the queue (using ``qstat``, falling back to ``qstat -H -f`` to
  determine whether a finished job succeeded or failed, since PBS's own
  job STATUS does not distinguish the two -- see Section 6)
- stop the workflow gracefully if a job fails and cannot be retried

**Run the master job from an interactive ``interact-g`` session, using
``--direct``.** This is the validated, recommended way to run SeisFlows on
Miyabi-G:

.. code:: bash

    # from the login node, get an interactive Miyabi-G allocation
    qsub -I -q interact-g -l select=1 -W group_list=<your project code>

    # once your interactive session starts, from your SeisFlows working dir
    seisflows submit --direct

``--direct`` runs the master job's control loop (file management, task
dispatch, queue monitoring) directly in your interactive shell's Python
process, instead of submitting the master job itself as a separate
``qsub`` job. This matters on Miyabi because:

- The ``prepost`` queue (which might otherwise seem like the natural place
  to run a lightweight control process) does **not** provide the same
  environment as the Miyabi-G login node -- if you land on a Miyabi-C
  (``x86_64``) pre-post node while your Conda environment/binaries are
  built for Miyabi-G (``aarch64``), nothing will run (see Section 6 for
  what that failure looks like).
- ``interact-g`` is real Miyabi-G hardware, so it correctly inherits your
  aarch64 Conda environment interactively, without needing the
  ``conda_env``/wrapper-script machinery described in Section 2 at all
  (that machinery is still required for the *compute* jobs the master
  process dispatches via ``run()`` -- only the master job's own execution
  benefits from running directly).

.. warning::

    ``seisflows submit --direct`` blocks in your interactive shell for as
    long as the workflow takes to run (the master job's control loop
    currently polls the PBS queue synchronously and does not return until
    the whole workflow finishes). ``interact-g`` sessions are capped at
    **2 hours**, and the session ends if your SSH connection drops. For
    anything beyond a quick test, run it inside ``tmux``/``screen`` so it
    survives a dropped connection, and be aware that a workflow needing
    more than ~2 hours of wall-clock master-job time will outlive a single
    ``interact-g`` allocation. There is no fully automated, unattended
    long-running story for Miyabi yet -- if you hit this limit, you
    currently need to re-request an ``interact-g`` session and re-run
    ``seisflows submit`` (without ``--direct``, or with it once you have a
    new interactive session), relying on SeisFlows's existing checkpoint
    file (``sfstate.txt``) to skip already-completed steps rather than
    starting over. This is a known limitation, being tracked separately.


7. How this works (site-specific PBS behaviors)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Getting the ``Pbs``/``Miyabi`` system module working against real Miyabi
hardware surfaced several PBS behaviors required bespoke tweaks. 

None of this is required reading to *use* the module, but it is useful context 
if something breaks or you are extending ``system.Pbs``/``system.Miyabi``:

- **qsub does not reliably forward trailing command line arguments to the
  submitted script.** Unlike SLURM's ``sbatch``, some ``qsub`` builds scan
  the *entire* command line for dash-prefixed tokens rather than stopping
  at the first non-option argument, so anything appended after the script
  name (e.g. ``--workdir /path``) gets fed back into ``qsub``'s own option
  parser and rejected as an unrecognized option
  (``qsub: invalid option -- '-'``, spelling out the flag character by
  character). Every piece of information a submitted script needs
  (working directory, pickled function/kwarg file paths, etc.) is instead
  passed as real job environment variables via ``qsub -v``, and
  ``runscripts/pbs_entry_point`` reads them back out and translates them
  into the CLI arguments the underlying ``submit``/``run`` scripts expect.
- **qsub copies a submitted script into its own spool directory before
  running it** (e.g. ``/var/spool/pbs/mom_priv/jobs/``), so a running
  script cannot reliably locate its sibling scripts by self-referencing
  its own path (bash's ``dirname "${BASH_SOURCE[0]}"`` resolves to that
  spool directory, not the real ``runscripts/`` directory). The real path
  is instead passed explicitly as ``SEISFLOWS_RUNSCRIPTS_DIR`` via
  ``qsub -v``.
- **PBS jobs do not start in the submission directory**, unlike SLURM's
  ``sbatch``. ``runscripts/pbs_entry_point`` explicitly ``cd``s into
  ``SEISFLOWS_WORKDIR`` (also passed via ``-v``) before doing anything
  else, since SeisFlows falls back to ``os.getcwd()`` for a module's
  working directory if one is not explicitly given.
- **Array jobs must be submitted rerunnable.** ``qsub -J`` fails with
  ``cannot submit non-rerunable Array Job`` unless ``-r y`` is also given,
  if the site's default ``Rerunable`` job attribute is ``False`` (as it is
  on Miyabi). ``Pbs.run_call()`` always adds ``-r y`` for array
  submissions.
- **PBS's ``-J`` only accepts a single contiguous ``start-end[:step]``
  range**, not an arbitrary comma-separated list of indices the way
  SLURM's ``--array`` does (``qsub: illegal -J value``). This matters for
  the ``rerun`` retry mechanism, which resubmits only the specific task
  indices that failed -- ``Pbs.run()`` groups failed indices into maximal
  contiguous runs and submits one ``qsub`` call per run, rather than one
  call with a comma-separated list.
- **A finished job's STATUS does not indicate success or failure.**
  ``qstat``/``qstat -H`` report states like ``FINISH`` regardless of
  whether the job actually succeeded, so success/failure has to come from
  elsewhere. The Miyabi User's Guide documents ``tracejob`` for this, but
  in practice it proved unreliable for array subjobs specifically --
  ``tracejob`` could not find records for either an individual subjob ID
  or the bare array parent sequence number, even for jobs confirmed to
  have completed successfully. ``qstat -H -f <job_id>`` (a specific job
  ID is required -- plain ``qstat -f`` only sees currently active jobs)
  reliably returns a finished job's full attribute set, including a line
  of the form ``Exit_status = 0``, and is what ``Pbs.query_job_states()``
  actually uses.
- **``qstat -t``/``qstat -H -t`` list an array's "parent"/summary row**
  (e.g. ``123456[].opbs``, empty brackets) alongside the real (sub)job
  rows. That parent row is a container, not an actual executed job, so it
  has no ``Exit_status`` of its own -- ``Pbs`` filters it out, since
  otherwise it would always be misreported as failed.
- **Conda needs its shell hook initialized before ``conda activate``
  works**, even though Conda is loaded as a module by default on Miyabi.
  Batch jobs run in a non-interactive shell that never sourced
  ``conda init``, so ``conda activate`` fails with
  ``CondaError: Run 'conda init' before 'conda activate'`` unless the
  hook is initialized first (``eval "$(conda shell.bash hook)"``), which
  ``runscripts/conda_activate-miyabi`` does automatically.

8. Troubleshooting
~~~~~~~~~~~~~~~~~~~~

- **"PBS 'queue' must match one of [...]"** -- the ``queue`` parameter must
  be one of Miyabi's node-occupied use queues (see Section 4 above, or
  the ``Miyabi`` class docstring for the current list and their node/
  walltime limits).
- **"qsub: invalid option -- '-'" (or similar, mentioning single
  characters)** -- this means something is appending extra command line
  arguments after the submitted script name; PBS's ``qsub`` does not
  reliably forward those to the script (Section 6). If you see this,
  you're likely running a modified/older version of the ``Pbs``/``Miyabi``
  classes -- update to the current version, which passes all such
  information via qsub's ``-v`` option instead.
- **"qsub: cannot submit non-rerunable Array Job"** -- your site defaults
  the ``Rerunable`` job attribute to false; ``Pbs.run_call()`` should
  already add ``-r y`` for every array submission (Section 6) -- if you
  still see this, you're likely running an older version of the module.
- **"qsub: illegal -J value"** -- PBS's ``-J`` only accepts a single
  contiguous range, not a comma-separated list of indices (Section 6).
  This should only be possible today if ``rerun`` is set and you are
  running an older version of the module whose retry logic predates the
  contiguous-range fix.
- **A job clearly ran and produced correct output, but SeisFlows reports
  it as FAILED** -- if this happens on a *current* version of the module,
  check the individual job's log file in ``logs/`` for the real error
  first, then double check ``tasktime`` (see the warning in Section 4) --
  a job killed for exceeding its walltime looks identical, in the logs,
  to any other failure. If neither explains it, compare
  ``qstat -H -f <job_id>`` (the mechanism ``Pbs`` actually relies on, see
  Section 6) against what SeisFlows reported.
- **"CondaError: Run 'conda init' before 'conda activate'"** -- the batch
  job's shell never sourced Conda's shell hook. This is handled
  automatically by ``runscripts/conda_activate-miyabi`` (via
  ``eval "$(conda shell.bash hook)"``); if you still see this, confirm
  ``conda_env`` is set (Section 2 above) so that wrapper is actually being
  used, and that plain ``conda`` (not just ``python``) is on ``PATH`` once
  its module is loaded.
- **A batch job fails immediately with a Python error like "parameter
  file does not exist" or an import/module-not-found error, despite
  everything above being correctly configured** -- this is the signature
  of the job running in the wrong directory or environment; confirm you
  are on the login node/subsystem matching your build (Section 0) and
  that ``conda_env`` (Section 2) points at a real, activatable
  environment. If you are developing/modifying ``Pbs``/``Miyabi``
  yourself, see the ``SEISFLOWS_WORKDIR``/``SEISFLOWS_RUNSCRIPTS_DIR``
  discussion in Section 6.
- **``conda``/Python "wrong architecture" errors on a pre-post or login
  session** (e.g. a Python script's shebang failing to exec, or bash
  trying to interpret Python source as a shell script) -- you are almost
  certainly on a Miyabi-C (x86_64) node trying to run an aarch64-built
  Conda environment, or vice versa. Check your hostname (``miyabi-g*`` vs
  ``miyabi-c*``) against Section 0.
- **Jobs rejected at submission** -- check ``show_token`` for remaining
  allocation, and confirm your ``group`` project code is correct.
- **Multibyte character errors** -- see the note in Section 1; rename your
  working directory/job title to plain ASCII.
- For anything else, check the individual job log files in ``logs/``, and
  compare against ``qstat -H -t <job_id>`` / ``qstat -H -f <job_id>`` run
  directly on the login node.