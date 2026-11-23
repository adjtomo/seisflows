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
