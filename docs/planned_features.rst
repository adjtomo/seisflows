Planned Features
================

Features that have been thought through but are not implemented, with enough
background to pick them up again later. Each section says what is wanted, what
was found out while looking at it, a proposed design, and what is still to be
decided. The working notes of the sessions that these come from are kept in
``docs/claude_contexts/``.

.. contents::
   :local:
   :depth: 2

--------------

seisflows restore
-----------------

**Status:** designed, not implemented (November 2026). Decided to wait, nothing
in the code yet.

What is wanted
~~~~~~~~~~~~~~

``seisflows clean`` removes the large files of a working directory that can be
made again by running more simulations (see `seisflows clean
<command_line_tool.html#seisflows-clean>`__). A complementary ``seisflows
restore`` should bring a cleaned working directory back to a state from which
the workflow can carry on, so that a User can free the disk space of a
stopped workflow, and later pick it up where it was.

The first idea was a small workflow or script that restores the deleted files
by copying model files back. Looking at the code showed that this is only
true for the models: most of what ``clean`` removes can only be made by
running simulations, and the normal workflow already does that. So ``restore``
is mostly about making the **state file** agree with the disk.

What ``clean`` removes and records
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``seisflows clean`` (``seisflows/tools/clean.py``) removes the files in these
categories, and keeps the parameter file, state file, logs, ``output/``,
``scratch/optimize/`` (the optimization vectors), residuals, ``scratch/system/``
and every directory (so, the run directory of every event):

``forward_arrays``
    ``proc*_save_forward_arrays*``, ``proc*_absorb_*`` in the run directories,
    and the ``*_FWD_ARR`` directories of a noise inversion.
``vtk``, ``databases``
    ``proc*.vt?`` and the other ``proc??????_*`` files of ``DATA``,
    ``OUTPUT_FILES``, ``OUTPUT_FILES/DATABASES_MPI`` and the ``LOCAL_PATH`` of
    the Par_file.
``kernels``
    ``eval_grad/kernels``, ``eval_grad/gradient``, ``eval_grad/misfit_kernel``.
``models``
    ``eval_func/model`` and ``eval_grad/model``.
``traces``
    ``traces/syn``, ``traces/adj``, ``SEM`` and the synthetics in
    ``OUTPUT_FILES``.
``obs``, ``solver_logs`` (only if asked for)
    ``traces/obs``, and timestamps and other small solver output.

Every file that it removes is appended to ``scratch/system/clean_manifest.txt``
with its category, size and path (relative to the working directory), under a
``# cleaned <time>`` line. This is what ``restore`` has to know what was lost.

What the workflow already makes again
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These were found by reading the code, and are why ``restore`` can be small:

- **Run directories.** ``solver.setup()`` runs on every ``seisflows submit``
  (also when a workflow is resumed). Its ``_initialize_working_directories``
  makes the directory structure again and copies ``bin`` and ``DATA``.
- **Databases and forward arrays.** ``solver.forward_simulation`` runs the
  mesher (``xmeshfem``, ``xgenerate_databases``, see
  ``_fwd_simulation_executables``) every time, from the model that
  ``solver.import_model`` has just imported. So running the forward stage
  remakes the databases, forward arrays, synthetics and (with the
  preprocessing) the adjoint sources.
- **Observed data.** ``prepare_data_for_solver`` is in the run list of
  ``Forward.evaluate_initial_misfit``, so the observed data is linked or copied
  again by every iteration.
- **Models.** ``Inversion.evaluate_initial_misfit`` writes ``m_new`` to
  ``eval_grad/model`` before it simulates (when the iteration is not thrifty).
  ``Inversion.initialize_line_search`` and ``Inversion.update_line_search``
  write the trial model ``m_try`` to ``eval_func/model``. The vectors that these
  come from, ``m_new``, ``m_try`` and ``alpha``, are in ``scratch/optimize``,
  which ``clean`` does not touch.

So a workflow that was cleaned between iterations, or just after
``finalize_iteration``, does not need anything restored: ``seisflows submit``
runs the first stage of the next iteration, which makes what it needs. The one
exception is a thrifty inversion (below).

What breaks
~~~~~~~~~~~

A workflow that was stopped in the middle of an iteration has a state file
that says that a task is done, while the files that the next task needs are
gone. The state file (``sfstate.txt``) has a line ``<task name>: 1`` for a
finished task, ``0`` for a pending one and ``PENDING:<job id>:<attempt>`` for
one that is dispatched, and the lines ``iteration``, ``was_thrifty`` and
``iteration_start``. The line search keeps a key for each of its steps
(``evaluate_line_search_misfit_i<iteration>s<step>``). The tasks of an
iteration are in ``Inversion.task_list``.

.. list-table::
   :header-rows: 1
   :widths: 28 32 40

   * - First pending task
     - It needs
     - If it is gone
   * - ``run_adjoint_simulations``
     - Forward arrays, databases and adjoint sources in the run directories.
     - Run ``evaluate_initial_misfit`` again (set its state, and everything
       after it, to ``0``): it re-runs the forward simulations.
   * - ``postprocess_event_kernels``
     - ``eval_grad/kernels``
     - Run ``run_adjoint_simulations`` again, which needs the row above.
   * - ``evaluate_gradient_from_kernels``
     - ``eval_grad/gradient`` (and the misfit kernel)
     - Run ``postprocess_event_kernels`` again, which needs the row above.
   * - ``perform_line_search``
     - ``eval_func/model``, the trial model
     - Write it again from ``scratch/optimize`` (``m_try``). No simulation.
   * - ``finalize_iteration``, or the next iteration, of a thrifty inversion
     - The forward arrays of the last line search step, which the next
       iteration uses instead of its own forward stage
     - Set ``was_thrifty`` to ``0``, and if ``thrifty`` is still on, say that
       the iteration being finalized cannot be thrifty.

The rows chain: the earliest task whose inputs are missing, back through what
that task needs, is where to rewind to.

Proposed design
~~~~~~~~~~~~~~~

``seisflows restore`` does not run simulations. It makes the working directory
consistent, and ``seisflows submit`` then does the work, as for any resumed
workflow:

1. Load the workflow as ``seisflows submit`` does (``import_seisflows``) and
   read the state file, the manifest and what is on disk.
2. Walk the first pending task back through what it needs, using the table
   above. Check the disk, not only the manifest, as the files may have been
   removed by something else or put back by hand.
3. Write the models that are missing from the optimization vectors (the only
   thing that is restored by copying), and rewind the state file as far as the
   walk says.
4. Print what it did and what it costs (e.g., "re-running ``evaluate_initial_misfit``,
   54 forward simulations"), and keep a copy of the state file before changing it.

``seisflows restore --dry_run`` prints this plan and changes nothing, like
``seisflows clean --dry_run``.

Suggested place in the code:

- ``Inversion.restore(dry_run=False)`` in ``seisflows/workflow/inversion.py``,
  that returns the plan so that it can be tested. ``Forward`` and ``Migration``
  only say that there is nothing to restore (they have no iterations to
  resume). ``NoiseInversion`` has its own state keys (one per force) and the
  ``*_FWD_ARR`` directories, and needs its own look.
- The command ``restore`` in ``seisflows/seisflows.py``, next to ``clean``,
  with ``-n/--dry_run`` and ``-f/--force``.
- Tests like the thrifty ones in ``seisflows/tests/test_workflow.py``: a
  working directory made by ``tests/test_clean.py`` (it has a file for every
  category), a state file with a task pending, and the plan that is expected for
  each row of the table. Nothing needs SPECFEM.

Still to decide
~~~~~~~~~~~~~~~

- **Rewind, or run?** The idea was a small workflow that restores the files.
  The proposal above only prepares the state and leaves the simulations to
  ``seisflows submit``. A variant, ``seisflows restore --run``, could submit a
  mini workflow that runs the forward stage straight away. Start with the
  rewind, which is what ``submit`` already does.
- **Thrifty.** If the forward arrays of a thrifty inversion were cleaned before
  ``finalize_iteration``, should ``restore`` set ``thrifty: False`` in the
  parameter file, or only say that it should be done? The suggestion is only to
  say it, as it should not edit parameters quietly.
- **A guard in ``clean``.** ``clean`` could warn when the state file shows an
  unfinished iteration, and name the categories (e.g., ``forward_arrays``) that
  would make ``restore`` rewind. It does not at the moment, the help only says
  that a workflow that carries on right after cleaning may need files that were
  removed.

What was not checked
~~~~~~~~~~~~~~~~~~~~

- The file patterns of ``clean`` come from the SPECFEM3D and SPECFEM3D_GLOBE
  code. They were **not** checked against a SPECFEM2D run directory
  (``DATA/proc*_*.bin`` as the models, and what SPECFEM2D puts in
  ``OUTPUT_FILES``), nor against a real working directory. Run ``seisflows clean
  --dry_run`` on one first.
- The table above says what each task reads, from reading the workflow code. It
  was not tested by cleaning a workflow stopped at each stage and resuming it.
- ``scratch/preprocess`` (Pyaflowa data and figures) is not touched by
  ``clean``, so ``restore`` does not consider it.

Related notes
~~~~~~~~~~~~~

``docs/claude_contexts/NONBLOCKING_MASTER_JOB_CHANGELOG.md`` has a section on
``seisflows clean`` with the same findings in more detail, and the thrifty
inversion is described in `Tips and Tricks <tips_and_tricks.html>`__.

--------------

Less time in the queue on busy clusters
---------------------------------------

**Status:** ideas, not implemented (November 2026). Written down after a
student's inversion on Miyabi took 171 hours for about 1 hour of compute.

Where the time goes
~~~~~~~~~~~~~~~~~~~

On a cluster with long queue times, the time of a workflow is the number of
trips through the queue, not the simulations. Every dispatch of a non-blocking
workflow waits twice: once as a compute job and once as the main job that
resumes the workflow. An iteration of an inversion is about 6 or 7 dispatches:
the initial misfit, the adjoint simulations, the kernel summing and smoothing
(a single-process job), and one dispatch for each step of the line search (the
unit step was rejected in 3 of the 4 iterations, so 3 or 4 of these).

Work that is done: ``pack``/``njob`` (see `cluster_setup
<cluster_setup.html#packing-tasks-into-fewer-jobs>`__) makes a stage wait once
instead of in waves of ``ntask_max`` jobs, ``partition_single`` sends the
single-process job to a short queue, and a thrifty inversion skips the initial
misfit from the third iteration on. The ideas below take trips away.

The three below were chosen to be written down. Others that were considered:

- **Better first line-search step.** Start at the last accepted step length (or
  a fraction of the unit step), so that the unit step is not rejected every
  iteration. Saves 1 or 2 dispatches an iteration, but it changes how the
  optimizer behaves, so it is a decision for the User.
- **Parallel line search.** Evaluate several trial steps (e.g., 1, 0.5, 0.25) in
  one dispatch and take the best that reduces the misfit. About 1.5 times the
  simulations of a sequential search but one trip instead of up to five. Needs
  the line search plugins to take several steps at once, and a solver directory
  for each step (so, not combined with thrifty).
- **A pilot job with a worker pool.** One long allocation runs a loop that picks
  up task files from the shared file system, so that after the first wait no
  stage waits in the queue. The idle time is only the serial steps between
  stages. The pack mode runner (``system/pack.py``) already runs each task as
  a subprocess with markers, which is most of a worker. The largest change and
  the largest win.

Merging the main job into the last compute job (so that there is no resume trip)
is a separate idea that is not written here.

Run the single-process stages in the main job
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**What is wanted.** Kernel summing and smoothing are dispatched as a
single-process job (``system.run(single=True)``), which is a trip through the
queue and a resume trip for about 4 minutes of work. The main job is already on
a compute node, so it could do the work itself and carry on, in blocking and
non-blocking mode alike. This saves about one dispatch (both trips) an
iteration.

**What was found.**

- The only callers with ``single=True`` are ``Migration.
  postprocess_event_kernels`` (``mask_source_event_kernels``,
  ``combine_event_kernels``, ``smooth_misfit_kernel``; inherited by
  ``Inversion``) and ``NoiseInversion``, which dispatches
  ``generate_event_kernels`` first (key ``postprocess_event_kernels_combine``)
  and then the one above.
- The functions are closures and methods of the workflow, so they can be called
  in the main process. ``Workstation.run`` already does this for a workstation:
  it sets task id 0 (``set_task_id``), runs the functions in this process and
  writes their output to a log file in ``logs/``, with a check of ``tasktime``
  after each function.
- ``Dispatcher.run`` is where a dispatch is chosen. A single-process stage that
  runs inline never suspends the workflow, so on non-blocking it just returns
  after the work, with no state key of its own to record (the task's key is
  written when it completes).
- They call SPECFEM binaries through ``solver.combine`` and ``solver.smooth``,
  which use ``system.mpiexec`` with ``nproc`` processes. A main job is
  submitted as one process: ``Pbs.submit_call_header`` has ``select=1:
  mpiprocs=1`` and ``Chinook.submit_call_header`` has ``--ntasks=1``. On Slurm
  ``srun`` cannot start more tasks than the job has. OpenMPI under PBS refuses
  more processes than slots. So the main job has to be asked for ``nproc``
  cores when this is on. Miyabi queues give whole nodes, so the cores are there,
  but the slots are not.

**Proposed design.**

1. A system parameter, e.g. ``single_inline`` (name to decide), default False,
   for ``Cluster`` systems. ``Dispatcher.run`` calls the ``Workstation.run``
   behaviour of the system for ``single=True`` when it is on, before any of the
   non-blocking code.
2. Main job size: when on, ``submit_call_header`` asks for ``nproc`` cores (PBS
   ``mpiprocs=min(nproc, node_size)``, Slurm ``--ntasks=nproc``), and the
   ``walltime`` of the main job has to cover the stage. Check that on
   ``check()`` (e.g., warn if ``walltime`` is less than ``tasktime``).
3. ``partition_single`` is not used when it is on. Say so in its docs.

**Cautions.**

- A main job that runs out of walltime in the middle of the stage ends without
  a resume job being queued (the next one is only queued by a dispatch), so the
  workflow just stops. ``seisflows submit`` carries on from the state file, but
  nothing says so. A warning when the stage starts, with the walltime that is
  left, would help.
- The smoothing can take much longer for large models (12 minutes for each
  parameter in the SPECFEM3D example with Gaussian smoothing), and a main job on
  ``debug-g`` has 30 minutes.
- Everything before the dispatch in ``postprocess_event_kernels`` runs again when
  the main job resumes, which is safe now and stays safe if the stage is
  inline (it removes the old ``misfit_kernel`` before it starts).

**Tests.** The simulated scheduler in ``tests/test_testflow.py`` and the
dummy workflows in ``tests/test_workflow.py`` can show that an inline stage
does not suspend, runs with task id 0, and does not dispatch. TestFlow could get
a check of its own (a single-process function that runs on the main job and
prints the number of cores it has), which is the part that needs a real
cluster.

Fuse the forward and adjoint simulations into one dispatch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**What is wanted.** ``evaluate_initial_misfit`` and ``run_adjoint_simulations``
are two dispatches, so two trips to the queue (and two resumes) an iteration.
Each event needs its own forward simulation before its adjoint one, and nothing
from the other events, so both can be one job.

**What was found.**

- The run list of ``Forward.evaluate_initial_misfit`` is ``prepare_data_for_
  solver``, ``run_forward_simulations`` and ``evaluate_objective_function``. The
  last one writes the adjoint sources of that event
  (``preprocess.quantify_misfit`` with ``save_adjsrcs``, in the run directory
  at ``traces/adj``), so after it the adjoint simulation of the event can run.
- ``run_forward_simulations`` saves the forward arrays for any workflow that is
  not ``Forward`` (``flag_save_forward``), and they stay in the run directory of
  the event, which is where the adjoint simulation reads them.
- ``Migration.run_adjoint_simulations`` dispatches ``_run_adjoint_simulation_
  single`` for each event, and ``Inversion`` overrides it to choose where the
  kernels are saved. The kernels summing is the first step that needs every
  event.
- A thrifty iteration skips the forward stage, so it has only the adjoint
  stage, as now.

**Proposed design.** A parameter of ``Migration``/``Inversion`` (e.g.,
``fuse_adjoint``, default False). When on, ``evaluate_initial_misfit`` appends
``_run_adjoint_simulation_single`` to its run list and ``run_adjoint_simulations``
is skipped (it checks its state key, or is removed from ``task_list``). The
``tasktime`` has to cover a forward and an adjoint simulation (about 3 + 5
minutes in the Miyabi inversion), which the docs have to say. A failed event
re-runs both. ``NoiseInversion`` has several forward simulations (a force each)
before its adjoint ones and is left out.

**Tests.** The run list that each dispatch gets, with the dummy workflows of
``tests/test_workflow.py``, with and without thrifty. Nothing needs SPECFEM.

Multi-node packed jobs
~~~~~~~~~~~~~~~~~~~~~~

**What is wanted.** On Miyabi the limit that binds is the number of jobs that
run at once for a project (3 or 4, see ``qstat --limit``), not the nodes. A job
may ask for several nodes: 8 in ``short-g`` and 16 in ``debug-g``. A packed job
with several nodes could run several tasks at the same time, one on each node,
so that ``njob`` jobs of 8 nodes run 8 times as many tasks at once. A stage of
54 events of 4 minutes would take minutes of compute rather than an hour, and
could fit into the 30 minutes of ``debug-g``, whose queue is much shorter.

**What was found.**

- A task is already allowed more than one node: ``Pbs.run_call`` asks for
  ``select={nodes}:mpiprocs=...`` with ``nodes`` from ``nproc`` and the node
  size. The packed job runs its tasks one after the other, in
  ``system.pack.run_chunk``, a loop of subprocesses.
- The solvers start their binaries with ``system.mpiexec`` (``solver.
  _run_binary``), which starts the processes on the nodes of the job, so tasks
  that run at the same time would all start on the same ones. Each needs its own
  nodes: a host list from ``$PBS_NODEFILE`` (or ``scontrol show hostnames
  $SLURM_JOB_NODELIST``) passed to ``mpirun``'s ``--host``/``-hostfile`` (or
  ``srun --nodelist``), by something the solver can read (e.g., an environment
  variable that ``mpiexec`` is built from).
- The GPU nodes of Miyabi have one GPU, so one task for each node.
- The markers and logs of ``pack`` do not change.

**Proposed design.** A parameter (e.g., ``nodes_per_job``, with ``njob``).
``run_chunk`` starts as many tasks at once as there are groups of ``nodes`` in
the job, each with its own hosts, and gives the next task of the block to the
first group that is free. The time limit of a job becomes ``tasktime`` times the
number of rounds (the block size divided by the groups, rounded up), plus the
buffer. Tasks in a block pick up work in turn, so the block sizes can stay even.

**Cautions.**

- Find out the node limit of the project first (``qstat --limit``). The
  Miyabi guide's table of limits was ambiguous about it (the number of nodes
  that the project can use at once, and per job).
- It uses more node-hours at once, which tokens are charged for, but the same
  number in all (tasks times node-hours of a task), apart from nodes that wait
  for the last task of a block.
- How ``mpirun`` of Miyabi behaves with a host list inside a PBS job has to be
  tried: this is the part that has to be tested on the cluster. TestFlow could
  print the hostname of each task, and check that tasks that run together are on
  different nodes.
