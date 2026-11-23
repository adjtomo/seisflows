# SeisFlows Testing Framework — Status

**Worktree**: `/Users/prof/Repos/seisflows/.claude/worktrees/testing-framework`
**Branch**: `feature/test-framework` (branched from `main`)
**Env**: micromamba env `seisflows` (already has the full runtime deps;
`pytest`, `pytest-cov`, `pytest-mock` were added to it during this work)

To open this worktree in VSCode: `code /path/to/.../testing-framework`
(don't try to check out this branch from the main repo's VSCode window —
it's already checked out here, and git will refuse).

## The ask

No robust way existed to test SeisFlows locally with confidence that a
change hasn't altered its *scientific* behavior. Two things make this hard:
SeisFlows drives real HPC job schedulers (SLURM/PBS/etc.), and most of its
actual compute work happens in external SPECFEM binaries. The `sfexamples`
(`seisflows/examples/*.py`) are the closest thing to an integration test
today, but they're demo scripts — not one of them asserts a misfit,
gradient, or model-update value. The goal: a layered testing framework that
gives real confidence, starting with pure unit tests (this session's focus)
and building toward scientific regression tests.

Proposed layers (cheapest/fastest first), for full context if resuming
this design conversation:
1. **Pure unit tests** — no subprocess, no SPECFEM. *(This session's work.)*
2. **Mocked-scheduler tests** — `unittest.mock.patch("subprocess.run")` to
   feed canned `sacct`/`qstat` output into `Slurm`/`Pbs`/`Miyabi`, so job
   query/retry/dependency logic is exercised without a real cluster.
   **Not started.**
3. **Local end-to-end workflow tests** — real `Forward`/`Inversion` task
   lists through `Workstation`/`Cluster` against a minimal solver, no
   scheduler. **Not started.**
4. **Golden-value regression tests** — a deliberately shrunk SPECFEM2D case
   run through real forward+adjoint+line-search, checked against stored
   reference misfit/gradient/model values with tolerance. This is the only
   layer that actually catches "the physics changed." **Not started.**

There is currently **no CI** in this repo (`.github/` has only an issue
template; `.travis.yml` is dead Python-2.7 cruft referencing a `scripts/`
layout that no longer exists — worth deleting whenever CI is wired up).

## What was done this session

### 1. Restored the existing suite to green
`seisflows/tests/*.py` predates the v3.5.5 "New Dev Strategy" rewrite
(commit `21e2b99b`) and had never been updated — it was **10 failed + 9
errored out of 35 tests** before this session touched anything. Every
failure was fixed against the *current* API (not worked around):

- `test_optimize.py` — `Model.save()` no longer exists (`Model.write()`
  only writes real SPECFEM Fortran-binary format now, no generic npz path);
  `Optimize.initialize_search()`/`update_line_search()` were split into
  separate primitives (`initialize_search()` + `calculate_step_length()` +
  `compute_trial_model()`, and `update_search()` + `calculate_step_length()`
  + `compute_trial_model()`); `step_len_max`/`step_len_min` safeguards now
  default to non-`None` and needed disabling for the toy Rosenbrock problem;
  three full-inversion tests (`test_inversion_optimization_problem_with_*`)
  are now real disk-I/O-per-iteration and were reduced/marked `slow`
  (see below) rather than run at their original 250-iteration budget.
- `test_tools.py` — `Model.save()`→`.write()`, plus dtype fix for
  `pytest.approx` on `dtype=object` arrays.
- `test_system.py` — `Cluster.run()`'s dummy test functions needed
  `**kwargs` (every real task-list function in the codebase accepts
  `**kwargs` because `Cluster.run()` always injects `verbose`/`log_level`).
- `test_preprocess.py` — `Default.read`/`.write`/`.initialize_adjoint_traces`
  moved from instance methods to module-level functions; `Pyaflowa` no
  longer pre-discovers `_station_codes`/`_source_names`/has `set_config()`
  in `setup()` — that responsibility moved to the calling workflow
  (`solver.source_names`), and `Pyaflowa.quantify_misfit()` now builds a
  per-call `Config` via `self._config.copy()` instead.
- `test_seisflows.py` — `SeisFlows.setup()` no longer creates the parameter
  file (`symlink`/`overwrite` kwargs are gone); that's now `init()`'s job.
  Renamed the stale test to `test_cmd_init` and added a lighter
  `test_cmd_setup_reports_parameter_errors` for current `setup()` behavior.

### 2. Two real bugs found along the way
- **`seisflows/optimize/NLCG.py` `checkpoint()` was completely broken** —
  it mutated a read-only `NpzFile` and referenced an undefined `dict_out`
  variable. Any real workflow using the NLCG optimizer would crash on
  `setup()`. **Fixed** (one-line, matches the working `LBFGS` pattern).
- **Environment gotcha, not a code bug**: the `seisflows` conda env's
  editable pip install points at `/Users/prof/Repos/seisflows` (the main
  checkout, currently on `feature-system_miyabi`), not this worktree. Plain
  `import seisflows` gets shadowed to the correct local copy by cwd-based
  `sys.path` insertion, but `Cluster.run()` spawns each task as a **separate
  subprocess** (via `system/runscripts/run`'s shebang) that does *not* get
  that shadowing — it was silently running the *other worktree's* code.
  Added `seisflows/tests/conftest.py` to fix `PYTHONPATH` for the whole
  pytest session. This will matter again for any Layer-2 `Slurm`/`Pbs`
  subprocess-dispatch tests, and for testing the `Pbs`/`Miyabi` work in
  general from a secondary worktree.

### 3. New coverage added (previously zero)
- **`seisflows/tests/test_workflow.py`** (new file, 12 tests) — covers
  `Forward`/`Inversion`'s state-file checkpoint/resume/failure control flow
  (skip completed tasks, mark failed tasks `-1` and stop, resume re-runs
  only failed/pending tasks, `stop_after`, per-iteration state reset,
  `is_thrifty`, `evaluation` string formatting). Uses dummy task lists via
  subclassing so no real `system`/`solver`/`preprocess`/`optimize` module is
  needed — this is the exact mechanism that would matter again if the
  non-blocking master-job redesign (see below) is picked back up.
- **`tools/noise.py` rotation math** (5 new tests in `test_tools.py`) —
  closed-form checks of the N/E↔R/T rotation matrices at simple angles.

### 4. Fast/slow tier split
Three Rosenbrock full-inversion tests and one real-Pyflex/Pyadjoint
`Pyaflowa` test take minutes each (real disk I/O per iteration / real
signal processing, not mocked). Marked `@pytest.mark.slow` and excluded by
default via `addopts = '-m "not slow"'` in `pyproject.toml`.

**Current state**: `pytest seisflows/tests` → 49 passed in ~2 min.
`pytest seisflows/tests -m slow` → 4 passed in ~10 min. **Nothing has been
committed yet** — everything above is uncommitted working-tree changes on
`feature/test-framework`.

## Where to resume

Pick up with whichever feels most valuable:
- **More Layer-1 coverage**: `Migration`/`NoiseInversion` still have no
  dedicated tests (only exercised indirectly via `Inversion`); `Slurm`'s
  pure-computation methods (`task_ids()`, `run_call()` string assembly,
  `_stdout_to_job_id()`) need zero mocking and are easy wins.
- **Layer 2**: mocked-scheduler tests for `Slurm`/`Pbs` using
  `unittest.mock.patch("subprocess.run")` with canned `sacct`/`qstat`
  output — this is the layer needed before touching the master-job
  non-blocking redesign with any confidence.
- **Layers 3/4**: the shrunk-SPECFEM2D golden-value regression harness —
  highest value for "did the physics change," highest effort to build.
- **Housekeeping**: decide whether/how to commit the current fixes,
  consider deleting the dead `.travis.yml`, and eventually wire up GitHub
  Actions (none exists today).

## Related, separate thread

This session started as a design discussion for making SeisFlows's master
job non-blocking on Miyabi (JCAHPC) — submit a job, then a scheduler-native
`afterany`-dependent job that re-invokes `seisflows submit`, so the master
process never needs to block/poll for a whole workflow's duration. That
design work is paused (not started in code) in favor of this testing
framework, which the user wants as a regression safety net *before* making
that kind of change. A separate, already-working `seisflows/system/pbs.py`
+ `seisflows/system/miyabi.py` module was being built independently in the
main checkout (`/Users/prof/Repos/seisflows`, branch
`feature-system_miyabi`) — check whether that has merged before resuming
either thread.
