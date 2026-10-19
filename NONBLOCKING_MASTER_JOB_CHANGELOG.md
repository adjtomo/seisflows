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

**Current total:** `pytest seisflows/tests` → 206 passed, 4 deselected (slow)

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
  - `WorkflowSuspended` (new exception, now in `dispatch.py`) — the
    sentinel raised after dispatching non-blockingly. Contract: persist
    `PENDING:<job_id>` as the task's state, checkpoint, call
    `self.system.submit(depends_on=job_id)`, then raise it. `Dispatcher.run`
    does all of this; tasks never raise it by hand.
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
- **State-file hardening** (prompted by a review of the checkpoint system
  before building Phase 1 on top of it — found that a single stray blank
  line in the file crashed the whole workflow load with a raw `ValueError`,
  confirmed by reproducing it):
  - `tools/config.py`: new `atomic_write_lines(path, lines)` — writes to a
    temp file in the same directory then `os.replace()`s it over the
    target, so a process killed mid-checkpoint (scheduler walltime limit,
    OOM) can't leave a truncated/unreadable state file. `Forward.checkpoint()`
    and `Inversion.checkpoint()` both switched to use it instead of
    overwriting in place.
  - `Forward.__init__`'s state-file reader now tolerates blank lines
    (harmless, e.g. a stray trailing newline from hand-editing) and raises
    a clear `logger.critical` + `sys.exit(-1)` naming the file and exact
    line for genuinely malformed content, instead of an unguarded
    `ValueError`/`int()` traceback.
  - Fixed the state file's self-written header comment, which still only
    documented `0`/`1`/`-1` and had gone stale the moment `PENDING:<job_id>`
    was introduced earlier in this same project.
  - Also fixed a latent one-character bug noticed in passing:
    `Inversion.checkpoint()`'s appended `iteration: N` line was missing its
    trailing `\n`.
  - 6 new tests: `test_tools.py` (`atomic_write_lines` round-trip, no
    leftover temp file, original file preserved if the replace step fails)
    and `test_workflow.py` (blank-line tolerance, clear error on a
    malformed value).
  - **Deliberately not done**: a lock/PID-stamp guard against two processes
    touching the state file concurrently. The non-blocking hand-off is
    sequential by design (dispatch → exit → *then* the dependent job can
    even be released), so this shouldn't occur in normal operation — the
    realistic risk is a user manually re-running `seisflows submit` while a
    dependent job is still pending. Revisit only if real-cluster testing on
    Chinook actually surfaces this.
- **Phase 1 pilot: non-blocking dispatch wired into real tasks.**
  (Originally written as methods on `Forward`; since moved into a
  `Dispatcher` class, see the "Dispatcher extracted" entry below. Names in
  this entry are the current ones.)
  - `Dispatcher.run(funcs, **kwargs)` (`seisflows/workflow/dispatch.py`) —
    the one call that task-list functions use instead of `system.run()`, as
    `self._dispatcher.run(...)`. Blocking systems: identical to
    `system.run()`. Non-blocking (`system.blocking == False`): the calling
    task is simply *entered twice*. First entry dispatches
    (`run(blocking=False)`), records `PENDING:<job_id>:<attempt>` as the
    task's state, checkpoints, schedules the dependent resume job, and
    raises `WorkflowSuspended`. Re-entry (in the resumed process) calls
    `system.job_outcome(job_id)`: success → return, so the code *after*
    the call (the "finish half") runs; failure → re-dispatch only the failed
    task ids and suspend again, up to `system.rerun` times, then mark the
    task `-1` and exit (so fixing the cause and re-running `seisflows
    submit` dispatches fresh, as after a blocking failure).
    **This closes the failure-handling gap** listed under the previous
    "Not done" — retries are honored across process restarts because the
    attempt count lives in the state file.
  - State value is now `PENDING:<job_id>:<attempt>` (attempt optional when
    parsing, so the earlier `PENDING:<job_id>` form still reads as attempt
    0). Parsed/validated by `parse_pending_state`; a malformed value gets
    the same clear state-file error as any other bad line.
  - Converted call sites: `Forward.generate_synthetic_data`,
    `Forward.evaluate_initial_misfit` (so `Inversion`/`NoiseInversion`
    overrides that call `super()` compose: on resume the whole override
    re-runs, and its post-processing after `super()` is the finish half),
    `Migration.run_adjoint_simulations`, and
    `Migration.postprocess_event_kernels`.
  - `Dispatcher.resuming`: True when the running task is re-entering to
    collect a dispatch. Needed because re-entry re-runs everything before
    the dispatch call: `postprocess_event_kernels` deletes the
    `misfit_kernel`/`mk_nosmooth` dirs *before* dispatching, which on resume
    would have destroyed the finished job's output, so that cleanup is now
    skipped when resuming. **Any task converted in future must be audited
    for this** (a pre-dispatch step that isn't safe to repeat).
  - Safety rails: at most ONE `dispatcher.run` call per task *and key* per
    invocation (a single PENDING state can't tell two dispatches apart) —
    a second raises `AssertionError` rather than silently skipping it on
    resume (tasks that dispatch several times give each its own `key`, see
    the Phase 2 entry). `_supports_nonblocking` (True on `Forward`, so on
    every shipped workflow) lets a workflow opt out, and `Forward.check()`
    then refuses `system.blocking=False` for it, instead of half-working.
    Any new workflow with a direct `system.run()` must set it False.
  - `system.blocking` (new `Cluster` param, default `True`) is the
    user-facing switch. `Cluster._supports_nonblocking` (False; True on
    `Slurm`) is declared per class, not inferred, and `Cluster.check()`
    refuses `blocking=False` where unsupported. Flipping the default to
    `False` is deliberately NOT done yet — wait for real-cluster validation.
  - `Slurm.job_outcome(job_id)` (generic `Cluster` stub raises
    `NotImplementedError`): reuses `monitor_job_status` (so it tolerates
    sacct lag) then reads failed task ids from each array element's own id
    (`'99_3'` → `3`) rather than list position — see Not done for why that
    matters.
  - **Resumed jobs append to the output log.** Found that `sbatch --output`
    truncates by default, and the Python logger writes to the same file
    (`sflog.txt`), so every resume would have erased the workflow's log
    history. Dependent submissions now add `--open-mode=append`
    (`Cluster._append_log_flag`, `Slurm` returns it) and skip the
    log/parameter-file backup `submit()` does for fresh submissions (it
    would otherwise pile up a numbered copy per resume). Systems with no
    append option fall back to backing up before each resume.
  - **A failed dependent submission is now fatal.** `Cluster.submit` ran
    `sbatch` without `check=True`, so a failed resume-job submission would
    have been silent: the master exits after dispatching and nothing ever
    resumes the workflow. It now exits with a message saying how to resume
    manually. (Fresh submissions are unchanged.)
  - Verified, not assumed: each resume re-runs `workflow.check()` and
    `workflow.setup()` (the `submit` runscript does) — checked every
    module's `setup()` and they are idempotent (mkdir / skip-if-exists;
    `Specfem._initialize_working_directories` skips existing source dirs),
    so a resume does not wipe a dispatched job's results.
  - Tests: `test_workflow.py` drives `Dispatcher.run` against a
    `_FakeSystem` (records dispatches/submits, replays scripted
    `job_outcome`s) — dispatch+suspend, resume success runs only the finish
    half, resume failure re-dispatches only failed ids, retry cap honored
    across restarts, fresh dispatch after exhaustion, second-dispatch
    assertion, `resuming` guard, `check()` gating. `test_system.py`:
    append/backup behavior, fatal failed resume submit, `job_outcome`
    (incl. partial re-run), support flags, `check()` gating.
- **Model validation moved off the master job (non-blocking mode only).**
  `solver.check_model_values` reads the whole model, which doesn't belong in
  a 1-core, short-walltime master (or on the login node under `--direct`).
  Split in two:
  - Master: `Specfem.check_model_files(path)` — file names only (path
    exists, every solver parameter has files, equal file counts across
    parameters). Never reads values; catches the common mistakes (wrong
    path, missing parameter, incomplete write) *before anything queues*.
  - Compute: `Specfem.validate_model_values(path)` — the expensive
    `Model.check()`; raises instead of exiting; stats printed by task 0
    only. Run as the first function of the dispatch
    (`Forward._check_model_single`, inserted at the top of the function list
    inline in each task, see the pattern below), told which model via
    `path_model`.
  - **Non-retryable failure.** A failed check writes its error to
    `<scratch>/preflight_failed.txt` and re-raises. After a failed job the
    master looks for that file *before* applying `rerun`: if present it
    exits with the check's own message, marks the task `-1`, and does NOT
    re-dispatch (a deterministic failure would just repeat and burn retries).
    The marker is cleared before every fresh dispatch and consumed when
    used, so a stale one can't be mistaken for a new failure. Scheduler
    agnostic (no exit-code parsing).
  - `Dispatcher.nonblocking` is the single definition of "dispatch is
    non-blocking right now" (a task is running AND `system.blocking` is
    False), shared by `Dispatcher.run` and the tasks' model checks.
  - The pattern, repeated inline in each task that checks a model (a few
    times, deliberately explicit rather than hidden behind a helper):
    ```python
    if self._dispatcher.nonblocking:
        self.solver.check_model_files(path=path_model)  # values on compute
        run_list.insert(0, self._check_model_single)
    else:
        self.solver.check_model_values(path=path_model)
    ```
  - Applied to `generate_synthetic_data` and `evaluate_initial_misfit`
    (`generate_synthetic_data` now passes `path_model` explicitly; same
    default it used before). **Blocking mode is untouched**: still the full
    `check_model_values` on the master, which stays the known-good baseline.
  - Costs to know about: in non-blocking mode every task reads the model to
    check it (they read it anyway in `import_model`, but watch the I/O at
    hundreds of tasks); and the model statistics now land in task 0's log
    under `logs/`, not in `sflog.txt`. A bad model that passes the file
    check now costs the queue wait before it is caught.
  - Tests: `test_workflow.py` (blocking vs non-blocking check placement,
    marker written/raised on failure, preflight failure not retried, other
    failures still retried, stale marker cleared), `test_solver.py`
    (structural check passes/catches missing path, missing parameter,
    incomplete write; never reads values; `validate_model_values` passes on
    a real model and raises rather than exits).
- **Decision: what NOT to move off the master.** Considered fusing
  master-side compute (e.g. gradient evaluation) into an adjacent
  `single=True` dispatch, and a separate preflight job with an `afterok`
  dependency. Both rejected: the master-side steps are small serial tasks, so
  a compute node is overkill and putting them on debug/prepost is too
  complicated; and `afterok` leaves the dependent array pending forever on
  failure unless the site sets `kill_invalid_depend`, with no clean PBS
  equivalent. Steps that mutate in-memory state (the optimizer) can't be
  offloaded transparently anyway. Optimizer and gradient evaluation stay on
  the master; if they prove too heavy, size the master through `walltime` and
  `slurm_args` (e.g. `--mem`) rather than moving them.

- **Dispatcher extracted from `Forward`.** The submission machinery was
  cluttering the workflow class, whose job is the science. It now lives in
  `seisflows/workflow/dispatch.py`:
  - `Dispatcher` (one per workflow, `self._dispatcher`, created in
    `Forward.__init__`): `run()`, `nonblocking`, `resuming`, the preflight
    marker (`marker`, `flag_preflight_failure`, `clear_preflight_failure`),
    and the suspend/resume/retry logic. `task(name)` is a context manager
    that `Forward.run()` wraps each task in, replacing the manual
    current-task/dispatch-count bookkeeping.
  - `WorkflowSuspended`, `pending_state()` and `parse_pending_state()` moved
    there too (`WorkflowSuspended` is still imported into `forward.py`).
  - `Forward` keeps `run()`, the `WorkflowSuspended` catch, the
    `_supports_nonblocking` gate in `check()`, and `_check_model_single`
    (it needs `self.solver`). `Migration` calls `self._dispatcher.run()` /
    `.resuming`.
  - The dispatcher holds a reference to its workflow and reads/writes task
    state through `workflow._states` / `workflow.checkpoint()` instead of
    owning any state. Chosen over passing getter/setter callables because
    `Inversion.run()` *reassigns* `self._states` each iteration (a held
    reference would go stale) and callables would add lambdas or tiny
    forwarding methods. The workflow<->dispatcher cycle pickles fine (tested:
    compute-side functions are bound methods of the workflow, which `dill`
    sends whole).
  - Removed `_with_model_check` / `_check_model`: the model-check branch is
    inline in `generate_synthetic_data` and `evaluate_initial_misfit`.
  - **The checkpoint system deliberately did NOT move into the dispatcher.**
    The state file records the progress of every task, blocking or not
    (`Forward.run()`'s skip-completed logic, `Inversion`'s iteration line,
    `NoiseInversion`'s sub-states); the dispatcher only needs to read and
    write one task's `PENDING` value. If it is ever extracted it should be
    its own object (a `Checkpoint`/state class owning read, validate,
    atomic write and reset) that both `Forward` and the `Dispatcher` take.
    Not done now: `_states` is touched in ~26 places across `forward.py`,
    `inversion.py` and `noise_inversion.py` plus the tests, which makes it a
    change of its own. The dispatcher's use of `workflow._states` is the one
    thing that would change.

- **Phase 2: the line search is a resumable loop; `Inversion` supports
  non-blocking.**
  - The recursive `update_line_search` could not be suspended: it applied
    `optimize.update_search()` (which consumes the saved `alpha`/`f_try`
    and advances `step_count`) BEFORE dispatching the next trial, so
    re-entering it after a resume would have applied that update a second
    time with the new `alpha` and the old `f_try`, corrupting the line
    search history.
  - New task `perform_line_search` replaces the `evaluate_line_search_misfit`
    and `update_line_search` tasks in `Inversion.task_list` (and
    `NoiseInversion.task_list`):
    ```python
    while True:
        self.evaluate_line_search_misfit()
        if self.update_line_search() == "PASS":
            break
    ```
    `update_line_search` is now only bookkeeping and returns `'PASS'`,
    `'TRY'` or `'RESTART'` instead of recursing; `evaluate_line_search_misfit`
    is one dispatch plus summing residuals into `f_try`. Both remain
    methods, only the task names changed.
  - Why resuming is correct: nothing is carried in memory between loop
    iterations. The optimizer checkpoint holds `step_count` and the search
    history (reloaded by `setup()` on every resume); `m_try`, `alpha` and
    `f_try` are vectors on disk. A resumed process re-enters at the current
    step: collect a pending evaluation, or dispatch the next. Every crash
    window is re-runnable: dying before the post-update checkpoint reloads
    the pre-update optimizer, so the update is recomputed; dying after it
    just dispatches the next step. `m_try` is written to `eval_func/model`
    before that checkpoint, so a step count in the checkpoint implies the
    model exists.
  - `Dispatcher.run(funcs, key=...)`: each dispatch is recorded under its
    own state-file key (default: the task name). The line search uses
    `evaluate_line_search_misfit_<evaluation>` (e.g. `..._i01s03`), so each
    step is tracked separately. Keys must not contain `:`. A task resets a
    key it wants to reuse: a line search FAIL restarts from step 1, so
    `update_line_search` zeroes every `evaluate_line_search_misfit*` key
    first. Keys left PENDING after being collected are harmless; they are
    reset with everything else at the end of the iteration.
  - The `m_try` model check follows the same inline `dispatcher.nonblocking`
    pattern, once per step.
  - Blocking mode uses the same loop (one implementation). Its crash
    recovery got simpler: a failed `perform_line_search` re-runs and
    re-enters at the persisted `step_count`, replacing the old trick of
    resetting `evaluate_line_search_misfit`'s task state to 0.
  - `Inversion` no longer opts out, so it passes `check()` with
    `blocking: False`. (`NoiseInversion` opted out at this point and was
    converted afterwards, see the next entry.)
  - **User-visible:** `stop_after` values `evaluate_line_search_misfit` and
    `update_line_search` are no longer valid, use `perform_line_search`. An
    inversion that is mid-iteration across this change re-runs its line
    search from the start (old task keys in the state file are ignored).
    Recorded in the release `CHANGELOG.md` under v3.6.0 (unreleased).
  - Tests (`test_workflow.py`, with a `_FakeOptimizer` that replays a script
    of statuses and, like the real one, only gives a resumed process its
    `step_count`): loop-until-PASS, PASS on step 1, suspend/resume across
    three processes (second collects only its own step, third does not
    re-collect step 1), no repeated bookkeeping on resume, failed step uses
    the dispatcher's rerun, FAIL-restart (keys reset, step 1 re-evaluated),
    FAIL with no restart exits, dispatcher keys (several dispatches in one
    task; `key` not passed on to the system), `check()` accepts
    `Inversion` and refuses a workflow that opts out.

- **`NoiseInversion` converted; every shipped workflow now supports
  non-blocking.** Its tasks dispatch once per force (and
  `postprocess_event_kernels` twice), so each dispatch gets its own state
  key and is skipped on resume once finished:
  - `generate_synthetic_data` (Z, E, N): key `generate_synthetic_data_<f>`.
  - `evaluate_rt_misfit` (E, N): key `evaluate_rt_misfit_<f>`. It reaches
    its dispatch through `super().evaluate_initial_misfit()`, so it uses the
    new `Dispatcher.keyed(name)` context manager instead of threading a
    `key` argument through `Inversion`/`Forward`'s signatures.
  - `_run_rt_adjoint_simulations_separate` (E/N x R/T): key
    `run_rt_adjoint_simulations_<f><c>`, again via `keyed()`.
  - `postprocess_event_kernels`: the horizontal-kernel combine uses key
    `postprocess_event_kernels_combine`; the smoothing in `Migration` keeps
    the task's own key, so its `resuming` guard (skip the cleanup that would
    delete the job's output) still works.
  - `evaluate_line_search_misfit` override: per-force keys
    `evaluate_line_search_misfit_<f>` (reset to 0 at the end of each step as
    before). Model check follows the inline `dispatcher.nonblocking`
    pattern; values are only checked before the first force of a step (and,
    for `generate_synthetic_data`, before the first force), since the rest
    only run once that has passed.
  - **A trap that was hiding here:** the per-force flags were tested with
    `bool(self._states[key])` / `in self._states and bool(...)`. A pending
    dispatch is a non-empty string, which is truthy, so the loop would have
    SKIPPED it (never collecting it) and moved on to the next force. Now
    `== 1`. A test confirms it fails if the old check is put back.
  - The per-force flags also make blocking crash-recovery finer: a resumed
    `evaluate_rt_misfit` / `generate_synthetic_data` skips forces already
    finished instead of redoing all of them.
  - **`Inversion._was_thrifty` is now persisted** in the state file as a
    `was_thrifty: 0|1` line next to `iteration`, written from the attribute
    in `Inversion.checkpoint()` (not from `_states`, which `run()` zeroes at
    the end of every iteration). Needed because `NoiseInversion.
    sum_all_residuals` reads it after dispatches that restart the process.
    See Not done for an apparent pre-existing problem with when it is reset.
  - `NoiseInversion` no longer opts out: `_supports_nonblocking` is True
    everywhere, with a test that guards it.
  - `workflow/test_flow.py` (`TestFlow`) was left alone at this point, and
    later converted, see "TestFlow can test a non-blocking main job".
  - Tests (`test_workflow.py`, real `NoiseInversion` methods against fakes,
    driving each task across several simulated processes): one suspend per
    force for synthetic data / RT misfit / line-search evaluation, the
    check placement (master: file names only, compute: first force only),
    the truthiness trap, the two keyed dispatches in `postprocess_event_
    kernels` including that its cleanup is skipped on resume, `keyed()`, and
    `was_thrifty` surviving a restart and the end-of-iteration reset.

### PBS / Miyabi port (`system/pbs.py`, written ahead of the Chinook run)

Done in parallel with, not after, the user's Chinook validation, since the
workflow side is scheduler-agnostic and this is the scheduler piece. It is
therefore the least-verified code here: mocked `qsub` only, and every
PBS-specific behaviour below is an assumption about the scheduler.

- `Pbs._supports_nonblocking = True` (inherited by `Miyabi`).
- `Pbs._dependency_flag` -> `-W depend=afterany:<id>`. `Pbs.submit(
  depends_on=)` appends it to the header; a failed resume submission is
  fatal and says how to resume by hand, as on SLURM.
- `Pbs.run(blocking=False)` returns the job id without monitoring. The
  qsub call moved into `_submit_run_call` (blocking path otherwise
  unchanged, still the baseline).
- **Retries need several jobs on PBS.** `-J` takes one contiguous range
  (SLURM's `--array` takes a list), so re-running failed tasks `1,3,4`
  submits `1-1` and `3-4`. `run(blocking=False, array="1,3,4")` therefore
  submits one job per contiguous range and returns their ids joined by
  `,`. That string is the dispatch's single job id as far as the
  `Dispatcher` and the state file are concerned (no change there: `,` is
  not the `:` the state file splits on), and `_dependency_flag` turns it
  into PBS's `a:b` for the resume job. If the second of two qsubs fails,
  the first is already running; it is fatal, resume by hand.
- `Pbs.job_outcome` monitors all the ids together and reads failed task ids
  from the subjob ids (`55[3].opbs` -> 3); a non-array job is task 0.
- `Cluster.monitor_job_status` list branch crashed (`+=` of `None`) when
  the scheduler had no answer yet for one job of a list. It now treats
  that as "no response yet" and keeps waiting. Nothing used the list form
  before this.
- **Logs differ from SLURM.** qsub has no append option, so there is no
  `_append_log_flag`; `Pbs.submit` already backed up the log before every
  submission, so each earlier master job's log is kept as a numbered copy
  in `path.log_files` and the live `sflog.txt` is replaced by the latest.
- Tests (`test_system.py`, mocked `subprocess.run`, on `Miyabi`):
  dependency flag incl. the `,` -> `:` conversion, `submit(depends_on)`
  with and without, fatal resume failure, `run(blocking=False)` incl.
  contiguous-range splitting and the single-task (no `-J`) case, blocking
  path unchanged, `_array_indices`, `job_outcome` for one job / several
  jobs / non-array, the monitor fix, and `check()` accepting
  `blocking: False`.
- **Things only Miyabi can answer** (check these first): is `-W depend=
  afterany:<id>` accepted alongside the existing `-W group_list=` (two `-W`
  options; if not they must be merged into one); does `afterany` on an
  array parent id (`55[].opbs`) wait for every subjob; does `qstat -H`
  still know a job when the resume job starts (history must be enabled, as
  the blocking path already assumes); are single-index arrays (`-J 3-3`)
  accepted (the blocking retry path already relies on it); is a dependent
  job accepted by `debug-g`; how a failed dependency behaves when the
  upstream job has already left the queue.

### Logging cleanup, phase 1: a resumed workflow reads as one story

Prompted by the first real run (2D example on Chinook, `sflog.txt` of nine
processes): every resume repeated ~17 lines of startup, banners and a
"skipping" line per finished task, about half the file. The user wants
`sflog.txt` to be the high-level story, with `logs/` for the detail and a
hint in the story about where to look. Agreed plan: phase 1 (this), phase 2
the INFO/DEBUG contract and an audit module by module, phase 3 new DEBUG
content and job queue/run times (`sacct`/`qstat`, user wants these).
Decisions: resume startup hidden entirely, "skipping" lines dropped.

- `Dispatcher.begin_resume()` (called at the end of `Forward.__init__`):
  if the state file has a PENDING dispatch, install `ResumeFilter` on the
  `seisflows` logger, which drops everything below WARNING. That hides the
  submit banner, `check()`/`setup()`, the iteration banner, the banners of
  finished tasks and whatever the pending task logs before its dispatch,
  with no change to any of those call sites. Warnings and errors still show.
- `Dispatcher.end_resume()` removes it, and is called when `run()` gets back
  to a dispatch (logging `---- RESUMING ITERATION 01: 'perform_line_search'
  after job 723824 ----`), at the first fresh dispatch, when `Forward.run`
  finishes and when a task raises. So the log can't stay hidden.
- "has already been run, skipping" is now DEBUG. `completed N tasks` is no
  longer logged when suspending (it said "0 tasks").
- On suspend, one INFO line replaces three: `dispatched job X, job Y will
  resume the workflow when it finishes. Task logs are in logs`. Needs the
  resume job's id, so `Cluster.submit`/`Pbs.submit` now capture the
  scheduler's response when resuming (it was the raw "Submitted batch job"
  line in the log) and return the id, logging the response at DEBUG.
- A failed task being re-run is a WARNING, naming the job and task ids.
- Demoted to DEBUG: "dispatched job X without waiting", "monitoring job
  status", "resuming key, checking outcome".
- Hints for where the detail is: `setup()` says where task output is
  written, the end of an iteration says where it was moved (`LOGS_iNN`,
  previously DEBUG), the suspend line names the log directory.
- Tests: the log of a resumed workflow is hidden until the dispatch is
  collected and then shown (fails with `begin_resume` disabled); a fresh one
  hides nothing; the filter is removed if the dispatch is never reached;
  the suspend line; `submit` returns the id and captures output only when
  resuming (Slurm and PBS). `conftest.py` clears logger filters after each
  test so a test that stops early can't hide another's log.
- Not yet: the `logs/` hint in failure messages, "running functions
  [...]" (INFO, unhelpful, wants phase 2), queue/run times.

### Logging cleanup, phase 2: the level contract, applied to the master job

The contract (INFO = the story, DEBUG = how, WARNING = carried on from,
CRITICAL = stopping) is written up in `docs/working_directory.rst`, under
"Which log to read". Applied here to what the master job writes to
`sflog.txt`, judged against the log of the first Chinook run. Message
levels only, plus the failure hints; no result values were added yet.

- INFO -> DEBUG, because they narrate a step whose result is logged next, or
  precede a check that logs nothing when it passes: `checking (true/target)
  model` (it logged a dangling colon), `m_try` check, `exposing gradient`,
  `calculating step length` / `computing line search trial model` /
  `calculating gradient dot products` (each followed by its result),
  `increment step count`, `saving misfit and step length` (the table after
  it stays INFO), the `?_new -> ?_old` / `m_try -> m_new` / `resetting step
  count` bookkeeping, `running finalization for module`, `cleaning out
  scratch directory`, `setting current iteration`, and `running functions
  [...]` in `Slurm`/`Pbs` (internal function names).
- Failure messages (`Dispatcher` retries spent, preflight failure, and the
  blocking `Slurm.run`/`Pbs.run`) now say `TASK LOGS: <dir>`.
- Left alone on purpose: compute-side messages (per source: `evaluating
  objective function for source`, ...), which are the story of `logs/`;
  the other systems' `running functions` (Cluster, Fujitsu, ...), which are
  blocking-only and the only record of what was run; everything in
  `optimize/` and `plugins/line_search/` apart from the above, which is
  already results and decisions.
- Test: a failure that stops the workflow names the task log directory.
- Next (phase 3, and what phase 2 leaves): add values where INFO announces a
  stage without a result (gradient min/max, a one-line iteration summary:
  misfit before/after, steps taken, wall time), job queue and run times,
  and an audit of `solver/`, `preprocess/` and `optimize/` DEBUG content,
  which at 51 `logger.debug` calls for the package is thin.

### Logging cleanup, phase 3: queue/run times and iteration summaries

- `Cluster.job_times(job_id)` (default None) -> `(queued, ran)` seconds.
  `Slurm`: `sacct -nLXP -o submit,start,end`, queue time is submit -> first
  start, run time first start -> last end over the array. `Pbs`: `qstat -H -f
  -t` per job (a re-run can take several), `qtime`/`stime`/`obittime` of each
  task. **The PBS version is a guess about what a PBS server reports** (the
  array itself has no times, `obittime` may be absent); it returns None
  rather than raising, so the worst case is the old log line. Check it on
  Miyabi. Both log the reason at DEBUG when they have nothing.
- `Dispatcher` logs them when it collects a job: `job 723809 finished
  successfully (queued 3s, ran 5s)`, and on a failure in the rerun warning.
  This is the number that answers whether resubmission latency matters.
- `msg.duration()` formats seconds (`3s`, `2m05s`, `1h02m`).
- `Inversion` logs `iteration 01 took 12m31s` after `COMPLETE ITERATION`.
  The start time is persisted as `iteration_start:` in the state file (like
  `was_thrifty`, since the process is restarted by a non-blocking
  dispatch), set when the iteration first runs and zeroed when it
  completes. The time includes queue waits, which is the point.
- `Gradient.finalize_search`: `misfit of accepted trial model is
  f=2.4E+01, 3.2% lower than before the line search (f=2.5E+01)`. It is
  the first element of the sorted search history, i.e. step length 0.
- Tests: `msg.duration`; Slurm and PBS time parsing incl. nothing to
  report; the dispatcher's line with and without times; the iteration start
  surviving a resume and the log line; the optimizer message.
- Not done: gradient min/max (the model update already logs `dm_min/dm_max`,
  so I held off), and the audit of DEBUG content in `solver/`,
  `preprocess/` and `optimize/`, which is the rest of the thin 51 calls.

### Logging cleanup, phase 3b: DEBUG content for the solver and preprocessing

These run in the tasks, so they end up in `logs/`, which is where users are
sent for detail. Judged by what is missing when a task goes wrong.

- `Specfem._run_binary`: DEBUG says the working directory and log file, and
  how long each executable took. **On failure the message now includes the
  log's real path (it named `scratch/solver/mainsolver/` whatever the
  source) and the last five lines of the solver's own log.** Before, finding
  out why SPECFEM failed meant going to look.
- `forward_simulation` / `adjoint_simulation`: DEBUG at the start (source,
  directory, `SAVE_FORWARD`), trace counts when exporting or moving, and the
  `FINISH ... SIMULATION` INFO says how long it took.
- `import_model`: DEBUG with the file count, and a **WARNING when it finds
  no model files** (previously a silent no-op that fails much later).
  `combine` says what it combines.
- `Default.quantify_misfit`: DEBUG with the number of obs-syn pairs and the
  components, a WARNING if there are none, `FINISH QUANTIFY MISFIT` now gives
  trace count, summed misfit and duration, and the per-trace `PREPROCESSING`
  line is DEBUG (one INFO line per trace was most of the file). A failing
  trace is named in `PREPROC FAILED`, with the traceback at DEBUG.
- Tests: the solver failure shows the log's path and last lines; importing
  nothing warns.
- Still to audit: `optimize/` DEBUG, `Pyaflowa`, the SPECFEM3D solvers'
  overrides of `forward_simulation`/`adjoint_simulation` (these changes are
  in the 2D/base class only), and the other systems.

### Detour: Example 4, a SPECFEM3D homogeneous halfspace (`seisflows examples run 4`)

Not part of the non-blocking work, but it is the natural larger-than-2D test
bed for it (`--par blocking=False`), so it lives here. Requested by the user
(2026-10-16): same approach as the 2D examples, a compiled SPECFEM3D is
pointed at (nothing is downloaded or built), initial model = the SPECFEM3D
`homogeneous_halfspace` example, target = the same with velocities +-5%, one
node (28 cores), `NSTEP` 2000, `XYZ` components, `--ntask` sources and
`--nsta` stations on a surface grid (default 25).

- `seisflows/examples/ex4_hh_specfem3d.py`, class `SFExample3D`, standalone
  (the 2D base class insists on a SPECFEM2D repo). Registered as choice 4 in
  `SeisFlows.examples()`; the listing finds it by file name.
- **Sizing:** `--nproc` (default 28 = a Chinook node) is split into the most
  square `NPROC_XI x NPROC_ETA` (28 -> 7x4) and `NEX` follows so elements stay
  about 3.7 km (28 -> 35x36). The mesh is made by `xmeshfem3D`, which is
  what the `Specfem3D` solver runs (the example's `run_this_example.sh` uses
  `xdecompose_mesh` instead, which SeisFlows does not).
- **Models:** `xmeshfem3D` + `xgenerate_databases` (MPI) write
  `proc*_vp/vs/rho/x/y/z.bin` because the example's Par_file has
  `SAVE_MESH_FILES`; they are copied to `OUTPUT_FILES_INIT`. The target is the
  same with the `Mesh_Par_file` material's vp and vs scaled (`--perturb`,
  default 5, density kept), then the original material is put back.
- **Sources / stations:** `CMTSOLUTION_001..N` (first at the center, the rest
  on a ring of 0.2 x the domain width), `STATIONS` a `ceil(sqrt(nsta))` grid
  with a 15% margin from the absorbing edges, equally spaced.
- **Cluster configuration:** `--system NAME`, `--par KEY=VALUE ...`,
  `--par_yaml FILE`. Parameters are validated against the generated parameter
  file first (`seisflows par` silently ignores an unknown key). Checked with
  `--system chinook --par partition=t1small submit_to=debug walltime=30
  tasktime=40 blocking=False`: all applied, `seisflows check` passes. Other
  new flags: `--specfem3d_repo`, `--perturb`.
- Refuses to run in a directory with an earlier workflow (`sfstate.txt`,
  `scratch`, `output`): a new parameter file would otherwise pick it up (the
  2D examples do, a stale state silently resumes).
- **Run for real, locally** (a throwaway SPECFEM3D build in the scratchpad,
  gfortran + Open MPI, the user's checkout untouched): `setup` and a full
  `run` with `--nproc 4 --nsta 9 --ntask 2`, one iteration, three line search
  steps, misfit down 94%. A 4-core forward simulation of 2000 steps took 32 s,
  the adjoint 1m30s. Also found by running it:
  - `smooth_type: gaussian` (the default) took **12 min per parameter** at 4
    cores; the example now uses `pde` (`xsmooth_sem_pde`, required binary).
  - the default `tasktime` of 1 minute timed the adjoint simulation out, the
    example now sets 30.
  - `Specfem3D.forward_simulation` wrote `ATTENUATION = .false`.` (a stray
    backtick); gfortran tolerated it. Fixed.
  - model generation needs a `DATA/CMTSOLUTION` and a pre-made
    `OUTPUT_FILES/DATABASES_MPI` (SeisFlows makes both for itself later).
- **Not validated:** a cluster. Nothing here has run under `srun`/`sbatch`,
  with 28 cores or with `blocking=False`; the next real run on Chinook is
  that test (from an interactive node: model generation is MPI). The
  default `walltime` of a cluster system (10 min) is for the non-blocking
  master job only; blocking needs a walltime for the whole workflow.
- Tests: `tests/test_examples.py` (28), with a stand-in SPECFEM3D repository,
  no SPECFEM needed: slices and `NEX` for several `nproc`, defaults, missing
  binaries / example, the working directory edits (Par_file, Mesh_Par_file,
  comments not mistaken for tables), ring sources, station grids, the
  target scaling, parameter validation and override order, the earlier
  workflow guard, the listing, and the argparse flags.
- Docs: `docs/specfem3d_example.rst`, in the examples toctree.

### TestFlow can test a non-blocking main job (no SPECFEM needed)

Requested by the user (2026-10-19): check the *system* module's non-blocking
parts on a cluster without compiling SPECFEM there.

- `TestFlow` is now a `Forward` subclass, so it gets the state file, the
  `Dispatcher` and the suspend/resume loop of every other workflow. Its
  tests call `self._dispatcher.run()` instead of `system.run()`, so the same
  five tests run blocking (as before) or with `blocking: False`.
  Set `seisflows par blocking False` (+ `submit_to`) and `seisflows submit`.
- Re-entry rules apply (see the dispatch docs): each test clears the log
  directory only `if not self._dispatcher.resuming`, and sets `rerun`
  inside the test (`setup()` runs again in every process and resets it).
  The checks of log files ignore `sflog*`/`parameters*` copies, which `Pbs`
  backs up into the log directory at every resume (qsub can't append).
- New test, non-blocking only: `test_preflight_failure_is_not_retried`. The
  task function writes the preflight marker to scratch and exits; the main
  job must stop without a re-run despite `rerun=2` and remove the marker.
  That checks the shared file system between compute and master as well.
- The hello world test also reads `system.job_times()` of its job (the job id
  comes from the state file): a warning, not a failure, if it is None. **This
  is the check for the guessed PBS fields.**
- "was expected to fail but did not" used to be a CRITICAL log message that
  carried on; it is an assertion now. Consequence: the generic local
  `cluster` system fails the single-job-failure test, because
  `Cluster.run(single=True)` ignores the return code of its one task
  (pre-existing, found by running it, not changed). `workstation` can't pass
  the array tests either (tasks run serially and stop at the first failure).
  The target is a real scheduler.
- When every test passes the state file is reset, so submitting again runs
  everything again (other workflows would skip finished tasks).
- **Bug found and fixed on the way: `dispatcher` was in every
  `parameters.yaml`.** `seisflows configure` writes each *public* attribute
  of a workflow as a parameter, and `self.dispatcher` (added in `70cfb1ae`)
  appeared as `dispatcher: <seisflows.workflow.dispatch.Dispatcher object at
  0x...>`. Harmless (the workflow ignores it) but wrong; the attribute is now
  `_dispatcher` everywhere and a test fails if a non-parameter becomes public.
  Existing parameter files can have the line deleted. The same mechanism made
  an early `self.system = ...` in `TestFlow.__init__` write `system: null`
  over the module choice, which is why it reads `self._modules` instead.
- Tests (`tests/test_testflow.py`): a `_SimulatedScheduler` stand-in that runs
  the dispatched tasks at once and remembers which failed. Blocking (one
  process, five jobs), non-blocking (seven processes, six jobs, a resume job
  waiting on each), the right tasks fail, it runs again after passing, the
  preflight test exists only non-blocking and runs one job despite `rerun`,
  the strict assertion, job times reported/missing, `__test__ = False`.
  Mutation checked: removing the `resuming` guard fails 7 of them.
- Not run on a real scheduler. First real use: `seisflows par workflow
  test_flow`, a cluster system, `blocking False`, on Chinook, then Miyabi.
  Expect ~6 trips through the queue (a few minutes each if it is busy).

## Not done / come back to this

- **Real-cluster validation — started, not finished.** First real run
  (2026-10-09, user, 2D example on Chinook, `t1small`, 1 source, one
  inversion iteration, `blocking: False`): **works.** Jobs queue behind the
  compute job with `afterany`, the workflow resumes exactly where it was,
  including inside the line search (three trial steps, each its own
  suspend/resume), the log appends across resumes, and the iteration
  completes. The user also reviewed the cleaned-up logging afterwards and
  was happy with it.
  Still to validate on a real scheduler: (1) a larger run (more `ntask`,
  larger simulations, real walltimes); (2) the failure paths, i.e. a
  deliberately failed array element to exercise `rerun` end to end, and a
  failed preflight (bad model) to check it is not retried; (3) the
  `forward`, `migration` and `NoiseInversion` workflows (only `Inversion`
  was run); (4) how big the per-resume overhead is at scale (startup was
  14-35 s per resume in the first run, about 8 s of it between the
  `workflow` and `system` setup lines, never profiled); (5) Miyabi/PBS,
  see below.
- **`Inversion._was_thrifty` was never reset — FIXED, see below.**
  `finalize_iteration` sets it True in its
  thrifty branch and nothing sets it back to False; its other branch wipes
  both `eval_grad` and `eval_func`. Sequence: iteration N is thrifty
  (flag True), N+1 is not (optimizer restarted in its line search, or it is
  the final iteration) so it is finalized the standard way and the flag
  stays True. If the inversion goes on to N+2 (a restart, not the final
  iteration), `evaluate_initial_misfit` sees the stale flag and skips the
  initial misfit, though `eval_grad` was just emptied. Fixed by setting
  `self._was_thrifty = False` in the standard `else` branch. This changes
  inversion behaviour independent of the non-blocking work, hence its own
  commit. `test_was_thrifty_follows_the_previous_iteration` runs the real
  `finalize_iteration` thrifty then standard, and failed before the fix.
- **Audit in-memory state that must survive a process restart: done for
  the shipped workflows.** A non-blocking workflow restarts the Python
  process at every dispatch, so anything a task stores on `self` and a
  later task reads is lost. Result: `Inversion._was_thrifty` was the one
  real case (now persisted, see the `NoiseInversion` entry);
  `NoiseInversion._force`/`_cmpnt` are set immediately before each use;
  `_srcrcv_stats` is rebuilt by `setup()` on resume; optimizer state is
  persisted via `optimize.checkpoint()`. **Re-check this for any task or
  workflow added later.**
- **PBS/Miyabi port — written, not validated on Miyabi** (see its entry
  above for the list of assumptions to check). Worth validating on Chinook
  first, since a workflow-side bug would otherwise be indistinguishable
  from a PBS-side one. The blocking retry path in `Pbs.run` has the same
  position-as-task-id flaw as `Slurm.run` (see below), left alone for the
  same reason.
- **`Frontera` gap — and it's already broken.** `Frontera(Slurm)`
  overrides `run()` itself, so it would not get `blocking=False` support for
  free (it declares `_supports_nonblocking = False` so `check()` refuses
  it). Separately, **`seisflows/system/frontera.py` can't be imported at
  all on this branch**: it does `from seisflows.system.slurm import Slurm,
  check_job_status`, and `check_job_status` doesn't exist in `slurm.py`
  (it's been that way since the v3.5.5 rewrite, `21e2b99b`). Pre-existing,
  unrelated to this project — left alone, but worth its own fix or removal.
- **Latent bug in the *blocking* retry path (not touched, on purpose).**
  `Slurm.run()` finds failed tasks by `enumerate()`-ing the sorted job list
  and treating the position as the task id. That's only right on the first
  attempt: after re-running just tasks `1,3`, the new array has two
  elements, so a repeat failure of task 3 is read as index 1 and the wrong
  task is re-run. The new `Slurm.job_outcome()` reads ids from the job ids
  instead. The blocking path is left untouched because it's the
  known-good baseline for validating the non-blocking one — once that's
  done, switch it to `job_outcome()` and delete the duplicated logic.
- **Master-job `walltime` guidance.** In non-blocking mode a master job
  only does bookkeeping, so `walltime` can be short (the point is to fit
  `debug`). Not yet documented for users; worth a docs note once validated.
  Measure first: the master still does gradient evaluation and optimizer
  vector operations, which scale with model size (see the decision above).
- **`Miyabi.submit_to`/`Chinook.submit_to`** (routing the master job to a
  cheap queue/partition) already exists and needs no changes — just a
  reminder that `submit_to="debug"` on Chinook is the real-world validation
  target once this works end-to-end, mirroring Miyabi's `debug-g` story.

## How to resume (if coming back cold)

1. Re-read this file in full before touching code.
2. `git log` on `feature/nonblocking-master-job` to confirm what's actually
   landed vs. what this file claims (this file could drift from reality).
   Last commit when this was written: `1c3b7a42` ("Add useful DEBUG
   content ..."). Nothing is pushed.
3. Where things stand: the non-blocking master job is built for every
   shipped workflow (`Forward`, `Migration`, `Inversion`, `NoiseInversion`)
   on SLURM/Chinook, and **validated on Chinook for a 2D `Inversion`**.
   PBS/Miyabi is written but only mocked. The logging cleanup (phases
   1-3b, see their entries) is done and the user has reviewed it.
4. **Next steps, in the order I'd take them.** The user drives, and is
   running more on Chinook, then will come back with instructions.
   - Whatever those runs turn up. Highest value: the failure paths and a
     larger run (see "Real-cluster validation").
   - Finish the DEBUG audit, started with `solver/specfem.py` and
     `preprocess/default.py`: `optimize/` DEBUG, `Pyaflowa`, the SPECFEM3D /
     3D_GLOBE overrides of `forward_simulation` / `adjoint_simulation`, the
     other systems. The user wants to say what they miss in DEBUG after
     real runs, so ask them first rather than guessing.
   - Logging odds and ends: gradient min/max in INFO was held back (the
     model update already logs `dm_min/dm_max`); revisit the contract in
     `docs/working_directory.rst` if INFO feels off.
   - Profile the per-resume startup (14-35 s) before it matters at scale.
     The `mkdir` loop over every path in `Forward.setup` is the first
     suspect, unconfirmed.
   - Miyabi/PBS validation (assumptions listed in the PBS entry, including
     the guessed `qstat` fields behind `Pbs.job_times`).
   - Once non-blocking is validated: switch the blocking `Slurm.run` /
     `Pbs.run` retry path to `job_outcome()` and delete the duplicated
     logic (latent position-as-task-id bug, see above), document master
     `walltime` guidance, and decide when `blocking` can default to False.
   - Parked, unrelated: `Frontera` doesn't import; optional `Checkpoint`
     class extraction (~26 `_states` sites).
5. User preferences (also in memory): pause before committing so they can
   review in VSCode, commit only when told; no single-line modifier helpers
   (inline, repetition is fine); keep blocking infrastructure out of
   `Forward`; 80-column max, PEP-8, docstrings with reST types; commit
   trailer `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`;
   never a bare `git stash`; remove stray `sflog.txt` before commits.
6. Tests: `micromamba run -n seisflows python -m pytest seisflows/tests`
   (206 passed, 4 slow deselected as of the last run).
