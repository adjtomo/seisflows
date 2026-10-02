# Non-Blocking Master Job — Changelog & Tracker

Internal tracker for the "make the SeisFlows master job non-blocking"
project. Update this file as work lands, not just at milestones, so a cold
resume (yours or someone else's) doesn't require reconstructing context from
commit messages or conversation history.

**Worktree**: `/Users/prof/Repos/seisflows/.claude/worktrees/testing-framework`
**Branch**: `feature/nonblocking-master-job` (renamed from
`feature/test-framework`, which this branch continues from — see
[TESTING_FRAMEWORK_STATUS.md](TESTING_FRAMEWORK_STATUS.md) for that
project's own history; it landed and is done, this file is for the new work)

## The problem

SeisFlows's master job blocks/polls (`Cluster.run()` → `monitor_job_status()`)
until each submitted job/array finishes, so it must live for the full
multi-hour-to-multi-day workflow. On Miyabi (JCAHPC) there's no acceptable
queue to host that cheaply (see the design discussion in conversation
history / ask the user for the original write-up if this needs
re-justifying) — `debug-g`/`debug-mig` (30 min max) is the target queue,
viable once the master job stops blocking.

## Design decisions (so they don't get re-litigated)

- **Scheduler-native dependency chaining, not polling removal.** Submit the
  compute job, then submit a second job with an `afterany`-style dependency
  on it, whose payload just re-invokes `seisflows submit`. The resume job's
  payload is the *existing* `submit()` mechanism, re-run — not a new entry
  point.
- **`afterany`, not `afterok`.** `afterok` would have SLURM auto-cancel the
  resume job the instant the upstream job fails (and for an array job,
  *any* element failing make the whole dependency unsatisfiable) — meaning
  the resumed process, which is what actually implements the existing
  `self.rerun`/retry-failed-array-indices logic, would never run. Failure
  handling has to stay SeisFlows's responsibility, not the scheduler's, so
  the dependency must hand control back unconditionally.
- **One implementation, not two maintained modes.** "Blocking" vs
  "non-blocking" is a thin branch inside `Cluster`/`Slurm.run()` (wait
  in-process vs. return `job_id` immediately), not two parallel workflow
  implementations. `blocking=True` stays the default everywhere (zero
  behavior change for existing callers) so it's nearly free to keep as a
  debug/validation baseline.
- **SLURM/Chinook first, PBS/Miyabi second.** User has better access to
  Chinook for real-world validation. The `workflow/*.py` dispatch/finish
  rewrite (the hard part) is 100% scheduler-agnostic and only needs doing
  once; the scheduler-specific piece (`_dependency_flag`, variable-passing)
  is thin and already structurally parallel between `Slurm`/`Chinook` and
  `Pbs`/`Miyabi` (see `pbs.py`'s existing `-v`-based variable passing vs.
  `Cluster.submit()`'s trailing-CLI-arg approach — each family reuses its
  own already-correct mechanism, just with a dependency flag added).
- **Rejected (for now): the "one big allocation, dispatch internally" idea**
  (user's proposal, matches long-standing GitHub issue #118 from 2022).
  Bigger/riskier than this project: requires writing a mini scheduler
  (bin-packing heterogeneous job-steps inside one allocation), risks
  *worsening* the original idle-node/resource-policy problem (not fixing
  it) unless paired with real cross-task backfill scheduling, has a much
  harder walltime-estimation problem (whole-iteration vs. per-task), and
  weaker PBS parity (`srun` job-steps have no clean PBS equivalent).
  Sequencing decision: finish non-blocking first — it'll tell us
  empirically whether resubmission latency is even a real problem worth
  solving the much harder way. Revisit only if needed.

## Done

- **`seisflows/system/cluster.py`**
  - `Cluster._dependency_flag(job_id)` — new, raises `NotImplementedError`
    (generic `Cluster` mimics a scheduler locally via `subprocess`/
    `ProcessPoolExecutor`; there's no real external job to depend on).
  - `Cluster.submit(..., depends_on=None)` — new param; when given, appends
    `self._dependency_flag(depends_on)` to the submit header. Ignored when
    `direct=True` (no scheduler call to attach a dependency to).
- **`seisflows/system/slurm.py`**
  - `Slurm._dependency_flag(job_id)` → `"--dependency=afterany:{job_id}"`.
    Inherited for free by `Chinook`, `Maui` (neither overrides `run()`/
    `submit()`).
  - `Slurm.run(..., blocking=True)` — new param, default `True` (no
    behavior change for existing callers). `blocking=False` submits and
    returns `job_id` immediately, skipping `monitor_job_status()` entirely.
    **Does not yet handle failure** — see Not Done.
- **`seisflows/workflow/forward.py`**
  - `WorkflowSuspended` (new exception) — the sentinel a task-list function
    raises after dispatching non-blockingly. Contract (see its docstring):
    persist `self._states[task_name] = f"PENDING:{job_id}"`, checkpoint,
    call `self.system.submit(depends_on=job_id)`, then raise it.
  - `Forward.run()` catches `WorkflowSuspended`, exits cleanly (not a task
    failure), and now **returns a status string** (`"COMPLETE"`/
    `"STOPPED"`/`"SUSPENDED"`) instead of `None`.
  - State-file reader fixed: was `line.split(":")` + `int(val)`, which
    breaks on a `PENDING:12345` value (extra colon, non-int). Now
    `split(":", 1)` with a conditional on the `PENDING` prefix.
- **`seisflows/workflow/inversion.py`**
  - `Inversion.run()` now checks `Forward.run()`'s returned status instead
    of the `self.stop_after is None` proxy it used before (which had a
    latent bug: if `stop_after` was set but never matched a task name, the
    outer loop would still stop after one iteration even though the task
    list ran to completion — fixed as a side effect). Does not advance
    `self.iteration`/reset states unless status is exactly `"COMPLETE"`.
- **Tests** (all passing, `pytest seisflows/tests` → 59 passed, 4
  deselected `slow`):
  - `test_system.py`: 7 new tests, mocked `subprocess.run` (no real
    cluster) — `_dependency_flag` on `Slurm`/`Cluster`, `submit(depends_on=
    ...)` building the right call (with/without dependency, `direct=True`
    ignoring it), `run(blocking=False/True)` returning `job_id`/`None` and
    calling (or not calling) `monitor_job_status`.
  - `test_workflow.py`: 3 new tests, dummy-subclass pattern (no real
    system/solver) — suspend/resume round-trip through the state file,
    `Inversion.run()` not advancing iteration on suspension.

## Not done / come back to this

- **Failure handling in `blocking=False` path — biggest gap.** Today,
  `run(blocking=False)` only covers the happy path. Nothing yet inspects
  whether the dispatched job actually succeeded. The resumed process (after
  the dependent `afterany` job fires) needs to replicate what
  `monitor_job_status()`'s blocking caller already does: query final job
  states, decide retry (respecting `self.rerun`/`_attempts`, resubmitting
  just the failed array indices) vs. give up with a clear error. Needs
  design: where does this check happen — generically in `Slurm` on resume,
  or per-task in each task-list function's "finish half"? Leaning toward
  the latter (consistent with the dispatch/finish split below), but not
  decided.
- **The actual dispatch/finish split in `workflow/*.py` — not started.**
  Everything above is plumbing; no real task-list function (e.g.
  `Forward.evaluate_initial_misfit`) has been rewritten to use it yet.
  Planned sequencing (from the original design discussion):
  - **Phase 1** (easier): non-recursive task-list functions in `Forward`,
    `Migration`, most of `NoiseInversion` — each calls `system.run()` once
    (or in a simple sequential loop) and does a deterministic thing with
    the result afterward. `NoiseInversion` already half-invented the
    needed idiom (per-substep `_states[key]` + `checkpoint()` in
    `_run_rt_adjoint_simulations_separate` and the `evaluate_line_search_
    misfit` override) for blocking crash-recovery — generalizing to the
    `PENDING` state is mechanical. Note: `NoiseInversion.
    postprocess_event_kernels` needs *two* sequential dispatches, not one.
  - **Phase 2** (hard): `Inversion.update_line_search()` is recursive — on
    "TRY" it calls `evaluate_line_search_misfit()` then itself,
    synchronously, deciding control flow from a result it just blocked to
    get. Must become an explicit iterative state machine (loop driven by
    state read at the top of each resume), not recursion.
  - Suggested pilot: `Forward.evaluate_initial_misfit` first (simplest
    real case), to validate the full contract end-to-end before touching
    anything else.
- **PBS/Miyabi port — not started, intentionally deferred.** Once validated
  on Chinook: add `Pbs._dependency_flag()` (`-W depend=afterany:{job_id}`)
  and a `blocking` param to `Pbs.run()`, reusing `Pbs.submit()`'s existing
  `-v`-based variable-passing (already correct for PBS's `qsub` quirks, see
  `pbs.py` module docstring) — same shape as the SLURM work, not a
  redesign.
- **`Frontera` gap.** `Frontera(Slurm)` overrides `run()` itself (own
  override, not inherited from `Slurm.run()`), so it does NOT get
  `blocking=False` support for free the way `Chinook`/`Maui` do. Needs its
  own matching update if/when non-blocking is wanted there. Not blocking
  for Chinook validation.
- **Real-cluster validation — not done.** Everything so far is mocked
  `subprocess.run()`; no real `sbatch --dependency=afterany` chain has been
  submitted to an actual SLURM scheduler. This needs to happen on Chinook
  before trusting the mechanism, even once Phase 1 lands.
- **`Miyabi.submit_to`/`Chinook.submit_to`** (routing the master job to a
  cheap queue/partition) already exists and needs no changes — just a
  reminder that `submit_to="debug"` on Chinook is the real-world validation
  target once this works end-to-end, mirroring Miyabi's `debug-g` story.

## How to resume (if coming back cold)

1. Re-read this file in full before touching code.
2. `git log` on `feature/nonblocking-master-job` to confirm what's actually
   landed vs. what this file claims (this file could drift from reality).
3. Check "Not done" above for the next concrete step — as of the last
   update here, that's the failure-handling design question, then the
   Phase 1 pilot.
