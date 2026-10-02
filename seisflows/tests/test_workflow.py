"""
Test the core workflow control-flow machinery in `workflow.forward.Forward`
and `workflow.inversion.Inversion` -- the state-file checkpoint/resume logic
that lets a workflow be safely re-submitted after a crash or job failure.

These tests deliberately avoid touching any real `system`/`solver`/
`preprocess`/`optimize` module: `task_list` is overridden with small dummy
functions so that the checkpoint/resume/failure bookkeeping in
`Forward.run()`/`Forward.checkpoint()` (and `Inversion`'s additional
iteration bookkeeping) can be exercised in isolation, quickly, with no I/O
beyond the state file itself.
"""
import os
import pytest

from seisflows.tools.config import Dict
from seisflows.workflow.forward import (Forward, WorkflowSuspended,
                                        _pending_state, _parse_pending_state)
from seisflows.workflow.inversion import Inversion


class _DummyForward(Forward):
    """
    Minimal `Forward` subclass with a 3-task list of no-op functions, so that
    `Forward.run()`'s checkpoint/resume/failure logic can be tested without
    any real system/solver modules.
    """
    def __init__(self, fail_on=None, suspend_on=None, **kwargs):
        self.calls = []
        self._fail_on = fail_on or set()
        self._suspend_on = suspend_on or set()
        super().__init__(**kwargs)

    @property
    def task_list(self):
        return [self.task_a, self.task_b, self.task_c]

    def task_a(self):
        self.calls.append("task_a")
        if "task_a" in self._fail_on:
            raise RuntimeError("task_a intentionally failed")

    def task_b(self):
        self.calls.append("task_b")
        if "task_b" in self._fail_on:
            raise RuntimeError("task_b intentionally failed")
        # Mimics a task-list function that calls `system.run(blocking=False)`,
        # gets back a job_id, and suspends the workflow to wait on it -- but
        # only the first time around; a resumed run (task_b already PENDING)
        # goes straight through, as the real `evaluate_initial_misfit`-style
        # dispatch/finish split would (finish-half logic isn't modeled here
        # since this harness only tests the suspend/resume control flow)
        if ("task_b" in self._suspend_on and
                self._states.get("task_b") != "PENDING:99999"):
            self._states["task_b"] = "PENDING:99999"
            self.checkpoint()
            raise WorkflowSuspended()

    def task_c(self):
        self.calls.append("task_c")


def test_forward_init_tolerates_blank_lines_in_state_file(tmpdir):
    """
    A stray blank line (e.g. a trailing newline added by a text editor) in
    the state file should be silently ignored, not crash the whole
    workflow -- it carries no data, unlike a genuinely malformed line
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")
        f.write("task_a: 1\n")
        f.write("\n")  # stray blank line

    workflow = _DummyForward(workdir=str(tmpdir), path_state_file=state_file)
    assert(workflow._states["task_a"] == 1)


def test_forward_init_raises_clear_error_on_malformed_state_file(tmpdir):
    """
    A non-blank line SeisFlows can't parse (e.g. a hand-edited value that
    isn't '0'/'1'/'-1'/'PENDING:<job_id>') should raise a clear, actionable
    error naming the file and line -- not a raw traceback from deep inside
    `int()`/tuple-unpacking
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")
        f.write("task_a: banana\n")

    with pytest.raises(SystemExit):
        _DummyForward(workdir=str(tmpdir), path_state_file=state_file)


def test_forward_run_executes_all_tasks_in_order(tmpdir):
    """A fresh workflow should run every task, in order, exactly once"""
    workflow = _DummyForward(workdir=str(tmpdir))
    workflow.path.state_file = os.path.join(str(tmpdir), "sfstate.txt")
    workflow.setup = lambda: None  # setup() needs real modules, skip it

    # `setup()` would normally create the state file; do the minimal
    # equivalent by hand
    with open(workflow.path.state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    status = workflow.run()

    assert(workflow.calls == ["task_a", "task_b", "task_c"])
    assert(workflow._states == {"task_a": 1, "task_b": 1, "task_c": 1})
    assert(status == "COMPLETE")


def test_forward_run_skips_completed_tasks_on_resume(tmpdir):
    """
    Re-running a workflow whose state file already shows every task as
    complete (`1`) should not re-execute any of them -- this is the
    mechanism that makes `seisflows submit` idempotent after a crash
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")
        f.write("task_a: 1\n")
        f.write("task_b: 1\n")
        f.write("task_c: 1\n")

    workflow = _DummyForward(workdir=str(tmpdir), path_state_file=state_file)
    workflow.run()

    assert(workflow.calls == [])


def test_forward_run_stops_and_marks_failed_task(tmpdir):
    """
    A task that raises should be marked `-1` (failed) in the state file, the
    exception should propagate (so the master job/scheduler sees a nonzero
    exit), and any tasks after the failed one should not run
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    workflow = _DummyForward(workdir=str(tmpdir), path_state_file=state_file,
                             fail_on={"task_b"})

    with pytest.raises(RuntimeError, match="task_b intentionally failed"):
        workflow.run()

    assert(workflow.calls == ["task_a", "task_b"])  # task_c never ran
    assert(workflow._states["task_a"] == 1)
    assert(workflow._states["task_b"] == -1)
    assert(workflow._states["task_c"] == 0)  # untouched, still pending


def test_forward_run_resumes_after_failure(tmpdir):
    """
    After a failed task is fixed and the workflow is re-submitted, already
    completed tasks should be skipped and the previously-failed task
    (state == -1, not 1) should re-run, along with anything after it
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    failing_run = _DummyForward(workdir=str(tmpdir), path_state_file=state_file,
                                fail_on={"task_b"})
    with pytest.raises(RuntimeError):
        failing_run.run()

    # Simulate a fresh process re-reading the same state file (as would
    # happen on a real `seisflows submit` resume), this time without the
    # induced failure
    resumed_run = _DummyForward(workdir=str(tmpdir), path_state_file=state_file)
    resumed_run.run()

    # task_a was already complete and should be skipped; task_b had failed
    # (-1, not 1) so it re-runs, and task_c runs for the first time
    assert(resumed_run.calls == ["task_b", "task_c"])
    assert(resumed_run._states == {"task_a": 1, "task_b": 1, "task_c": 1})


def test_forward_stop_after(tmpdir):
    """`stop_after` should end the workflow early without marking later
    tasks as failed -- they're simply never attempted"""
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    workflow = _DummyForward(workdir=str(tmpdir), path_state_file=state_file,
                             stop_after="task_b")
    status = workflow.run()

    assert(workflow.calls == ["task_a", "task_b"])
    assert(workflow._states["task_a"] == 1)
    assert(workflow._states["task_b"] == 1)
    assert(workflow._states["task_c"] == 0)
    assert(status == "STOPPED")


def test_forward_run_suspends_non_blocking_dispatch(tmpdir):
    """
    A task that raises `WorkflowSuspended` (mimicking `system.run(
    blocking=False)` dispatching work and the task persisting its own
    'PENDING:<job_id>' state) should stop the workflow cleanly -- not as a
    failure, and without touching tasks after it
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    workflow = _DummyForward(workdir=str(tmpdir), path_state_file=state_file,
                             suspend_on={"task_b"})
    status = workflow.run()

    assert(workflow.calls == ["task_a", "task_b"])  # task_c never ran
    assert(status == "SUSPENDED")
    assert(workflow._states["task_a"] == 1)
    assert(workflow._states["task_b"] == "PENDING:99999")  # not failed (-1)
    assert(workflow._states["task_c"] == 0)


def test_forward_run_resumes_pending_task_after_suspension(tmpdir):
    """
    Re-reading a state file with a 'PENDING:<job_id>' task should round-trip
    correctly (not raise/crash on `int()`-parsing a job ID), and a fresh
    process picking it back up should be able to tell 'this task already
    dispatched work' apart from 'this task hasn't started yet'
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    first_run = _DummyForward(workdir=str(tmpdir), path_state_file=state_file,
                              suspend_on={"task_b"})
    first_run.run()

    # Simulate the dependent job re-invoking `seisflows submit`: a fresh
    # process re-reads the same state file
    resumed = _DummyForward(workdir=str(tmpdir), path_state_file=state_file,
                            suspend_on={"task_b"})
    assert(resumed._states["task_a"] == 1)
    assert(resumed._states["task_b"] == "PENDING:99999")

    status = resumed.run()

    # task_a already complete (skipped); task_b's dummy "finish half" (see
    # `_DummyForward.task_b`) runs without re-suspending since it's already
    # PENDING; task_c runs for the first time
    assert(resumed.calls == ["task_b", "task_c"])
    assert(status == "COMPLETE")
    assert(resumed._states == {"task_a": 1, "task_b": 1, "task_c": 1})


class _DummyInversion(Inversion):
    """
    Minimal `Inversion` subclass with a single dummy task, used to test the
    per-iteration bookkeeping in `Inversion.run()`/`Inversion.checkpoint()`
    without needing real system/solver/preprocess/optimize modules
    """
    def __init__(self, suspend_on_iteration=None, **kwargs):
        self.calls = []
        self._suspend_on_iteration = suspend_on_iteration
        super().__init__(**kwargs)

    @property
    def task_list(self):
        return [self.only_task]

    def only_task(self):
        self.calls.append(self.iteration)
        if (self.iteration == self._suspend_on_iteration and
                self._states.get("only_task") != "PENDING:99999"):
            self._states["only_task"] = "PENDING:99999"
            self.checkpoint()
            raise WorkflowSuspended()


def test_inversion_run_increments_iteration_and_resets_states(tmpdir):
    """
    Each full pass through the task list should increment `self.iteration`
    and reset every task's state back to pending (0) for the next iteration
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    workflow = _DummyInversion(workdir=str(tmpdir), path_state_file=state_file,
                               start=1, end=3)
    workflow.iteration = workflow.start  # normally set by setup()
    workflow.run()

    assert(workflow.calls == [1, 2, 3])
    # After the final iteration, `run()` increments `iteration` past `end`
    # before the loop condition fails, and states are reset for a (nonexistent)
    # next iteration
    assert(workflow.iteration == 4)
    assert(workflow._states["only_task"] == 0)


def test_inversion_run_does_not_advance_iteration_on_suspension(tmpdir):
    """
    If a task-list function suspends the workflow mid-iteration (non-blocking
    dispatch), `Inversion.run()` must NOT increment `self.iteration` or reset
    states for a next iteration -- that would silently skip re-checking the
    in-progress (suspended) iteration when a dependent job resumes it
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    workflow = _DummyInversion(workdir=str(tmpdir), path_state_file=state_file,
                               start=1, end=3, suspend_on_iteration=2)
    workflow.iteration = workflow.start
    workflow.run()

    # Stopped mid-iteration 2, having completed iteration 1 in full
    assert(workflow.calls == [1, 2])
    assert(workflow.iteration == 2)  # NOT advanced to 3
    assert(workflow._states["only_task"] == "PENDING:99999")

    # A dependent job re-invoking `seisflows submit` resumes from here: a
    # fresh instance re-reads iteration==2 and the PENDING task state, picks
    # back up in the same (not-yet-advanced) iteration, and proceeds normally
    resumed = _DummyInversion(workdir=str(tmpdir), path_state_file=state_file,
                              start=1, end=3, suspend_on_iteration=2)
    assert(resumed.iteration == 2)
    resumed.run()

    assert(resumed.calls == [2, 3])
    assert(resumed.iteration == 4)


@pytest.mark.parametrize(
    "thrifty,iteration,start,end,restarted,expected", [
        (False, 2, 1, 5, False, False),  # thrifty disabled entirely
        (True, 1, 1, 5, False, False),   # first iteration never thrifty
        (True, 5, 1, 5, False, False),   # final iteration never thrifty
        (True, 3, 1, 5, True, False),    # optimizer restarted -> not thrifty
        (True, 3, 1, 5, False, True),    # otherwise: thrifty
    ])
def test_inversion_is_thrifty(tmpdir, thrifty, iteration, start, end,
                              restarted, expected):
    """
    `is_thrifty` should only be True when thrifty is enabled AND we're
    strictly between the first and last iteration AND the optimizer hasn't
    just been restarted
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    workflow = _DummyInversion(workdir=str(tmpdir), path_state_file=state_file,
                               start=start, end=end, thrifty=thrifty,
                               optimize="LBFGS")
    workflow.iteration = iteration

    class _StubOptimize:
        _restarted = restarted

    workflow.optimize = _StubOptimize()

    assert(workflow.is_thrifty == expected)


def test_inversion_evaluation_string(tmpdir):
    """
    `evaluation` should format `iteration`/`step_count` into the
    'i{iter}s{step}' tag used throughout the codebase for file naming
    """
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    workflow = _DummyInversion(workdir=str(tmpdir), path_state_file=state_file,
                               start=1, end=1)
    workflow.iteration = 3

    class _StubOptimize:
        step_count = 7

    workflow.optimize = _StubOptimize()

    assert(workflow.evaluation == "i03s07")


# -----------------------------------------------------------------------------
# `Forward._run_on_system()`: non-blocking dispatch, resume, and retry
# -----------------------------------------------------------------------------

class _FakeSystem:
    """
    Stand-in for a system module that records what the workflow asks of it
    and replays scripted job outcomes, so the dispatch/resume/retry logic can
    be exercised without a scheduler
    """
    def __init__(self, par_file, blocking=False, rerun=0, outcomes=None):
        self.blocking = blocking
        self.rerun = rerun
        self.path = Dict(par_file=par_file)
        self.outcomes = list(outcomes or [])  # [(success, failed_ids), ...]
        self.runs = []  # (blocking kwarg or None, array kwarg or None)
        self.submits = []  # depends_on of each resume job scheduled
        self._next_job_id = 1000

    def check(self):
        pass

    def run(self, funcs, blocking=None, array=None, **kwargs):
        self.runs.append((blocking, array))
        if blocking is False:
            self._next_job_id += 1
            return str(self._next_job_id)

    def submit(self, workdir=None, parameter_file=None, depends_on=None):
        self.submits.append(depends_on)

    def job_outcome(self, job_id):
        return self.outcomes.pop(0)


class _DummySystemWorkflow(Forward):
    """
    Workflow whose second task dispatches to the system through
    `_run_on_system()` and then does some 'finish' work, mirroring e.g.
    `Forward.evaluate_initial_misfit`
    """
    def __init__(self, **kwargs):
        self.calls = []
        super().__init__(**kwargs)

    @property
    def task_list(self):
        return [self.before, self.dispatching_task, self.after]

    def before(self):
        self.calls.append("before")

    def dispatching_task(self):
        self.calls.append("dispatching_task:start")
        self._run_on_system([self.before], x=1)
        self.calls.append("dispatching_task:finish")

    def after(self):
        self.calls.append("after")


def _make_workflow(tmpdir, **system_kwargs):
    """Build a `_DummySystemWorkflow` wired to a `_FakeSystem`"""
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    if not os.path.exists(state_file):
        with open(state_file, "w") as f:
            f.write("# SeisFlows State File\n")
    workflow = _DummySystemWorkflow(workdir=str(tmpdir),
                                    path_state_file=state_file)
    workflow.system = _FakeSystem(par_file=os.path.join(str(tmpdir), "p.yaml"),
                                  **system_kwargs)
    return workflow


def test_pending_state_round_trip():
    """PENDING values carry the job id and retry count, and parse back"""
    assert(_parse_pending_state(_pending_state("123", 2)) == ("123", 2))
    assert(_parse_pending_state("PENDING:123") == ("123", 0))
    assert(_parse_pending_state("PENDING:123[].opbs:1") == ("123[].opbs", 1))


@pytest.mark.parametrize("val", ["PENDING", "PENDING:", "PENDING:1:2:3",
                                 "PENDING:1:x", "OTHER:1:0"])
def test_parse_pending_state_rejects_malformed(val):
    with pytest.raises(ValueError):
        _parse_pending_state(val)


def test_state_file_with_malformed_pending_raises_clear_error(tmpdir):
    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("task_a: PENDING:\n")
    with pytest.raises(SystemExit):
        _DummyForward(workdir=str(tmpdir), path_state_file=state_file)


def test_run_on_system_blocking_is_plain_system_run(tmpdir):
    """With `system.blocking` True nothing about dispatch changes"""
    workflow = _make_workflow(tmpdir, blocking=True)
    status = workflow.run()

    assert(status == "COMPLETE")
    assert(workflow.calls == ["before", "dispatching_task:start",
                              "dispatching_task:finish", "after"])
    assert(workflow.system.runs == [(None, None)])  # no `blocking` kwarg
    assert(workflow.system.submits == [])


def test_run_on_system_outside_a_task_is_blocking(tmpdir):
    """
    Called directly (e.g. from `seisflows debug`) rather than from
    `Forward.run()`, there's no task to suspend, so it just blocks
    """
    workflow = _make_workflow(tmpdir, blocking=False)
    workflow._run_on_system([workflow.before])

    assert(workflow.system.runs == [(None, None)])
    assert(workflow.system.submits == [])


def test_run_on_system_nonblocking_dispatches_and_suspends(tmpdir):
    workflow = _make_workflow(tmpdir)
    status = workflow.run()

    assert(status == "SUSPENDED")
    # Stopped at the dispatch: nothing after it ran, not even the finish half
    assert(workflow.calls == ["before", "dispatching_task:start"])
    assert(workflow.system.runs == [(False, None)])
    # Resume job scheduled behind the dispatched job
    assert(workflow.system.submits == ["1001"])
    assert(workflow._states["dispatching_task"] == "PENDING:1001:0")
    # ...and the state file already says so, for the resumed process to read
    resumed = _make_workflow(tmpdir)
    assert(resumed._states["dispatching_task"] == "PENDING:1001:0")


def test_run_on_system_resume_success_runs_finish_half_only(tmpdir):
    _make_workflow(tmpdir).run()

    resumed = _make_workflow(tmpdir, outcomes=[(True, "")])
    status = resumed.run()

    assert(status == "COMPLETE")
    # `before` already complete (skipped); the dispatching task is re-entered
    # but does NOT dispatch again; its finish half and `after` then run
    assert(resumed.calls == ["dispatching_task:start",
                             "dispatching_task:finish", "after"])
    assert(resumed.system.runs == [])
    assert(resumed.system.submits == [])
    assert(resumed._states["dispatching_task"] == 1)


def test_run_on_system_resume_failure_reruns_only_failed_tasks(tmpdir):
    _make_workflow(tmpdir, rerun=2).run()

    # Resumed process finds tasks 1 and 3 failed: re-dispatch just those
    second = _make_workflow(tmpdir, rerun=2, outcomes=[(False, "1,3")])
    assert(second.run() == "SUSPENDED")
    assert(second.system.runs == [(False, "1,3")])
    assert(second.system.submits == ["1001"])
    assert(second._states["dispatching_task"] == "PENDING:1001:1")

    # ...and the retry succeeds
    third = _make_workflow(tmpdir, rerun=2, outcomes=[(True, "")])
    assert(third.run() == "COMPLETE")
    assert(third.system.runs == [])


def test_run_on_system_retries_are_capped_across_restarts(tmpdir):
    """
    The retry count lives in the state file, so `system.rerun` is honored
    even though every attempt happens in a different process. Once spent,
    the workflow exits and marks the task failed (not PENDING) so that
    fixing the cause and re-running `seisflows submit` dispatches fresh
    """
    _make_workflow(tmpdir, rerun=1).run()  # original dispatch
    second = _make_workflow(tmpdir, rerun=1, outcomes=[(False, "0")])
    assert(second.run() == "SUSPENDED")  # attempt 1 of 1

    third = _make_workflow(tmpdir, rerun=1, outcomes=[(False, "0")])
    with pytest.raises(SystemExit):
        third.run()
    assert(third._states["dispatching_task"] == -1)

    fresh = _make_workflow(tmpdir, rerun=1)
    assert(fresh._states["dispatching_task"] == -1)
    assert(fresh.run() == "SUSPENDED")
    assert(fresh.system.runs == [(False, None)])  # a fresh dispatch


def test_run_on_system_rejects_second_dispatch_in_one_task(tmpdir):
    """
    One PENDING state can't tell two dispatches apart, so a task making a
    second non-blocking dispatch must fail loudly rather than silently skip
    the second one on resume
    """
    class _TwoDispatches(_DummySystemWorkflow):
        def dispatching_task(self):
            self._run_on_system([self.before])
            self._run_on_system([self.before])

    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")
        f.write("dispatching_task: PENDING:1001:0\n")
    workflow = _TwoDispatches(workdir=str(tmpdir), path_state_file=state_file)
    workflow.system = _FakeSystem(par_file="p.yaml", outcomes=[(True, "")])

    with pytest.raises(AssertionError, match="more than once"):
        workflow.run()


def test_resuming_flag_guards_unsafe_preamble(tmpdir):
    """
    Work a task does before dispatching that destroys the dispatched job's
    output (e.g. `Migration.postprocess_event_kernels` clearing kernel
    directories) must be skipped when resuming
    """
    seen = []

    class _Guarded(_DummySystemWorkflow):
        def dispatching_task(self):
            seen.append(self._resuming)
            self._run_on_system([self.before])

    state_file = os.path.join(str(tmpdir), "sfstate.txt")
    with open(state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    first = _Guarded(workdir=str(tmpdir), path_state_file=state_file)
    first.system = _FakeSystem(par_file="p.yaml")
    first.run()

    resumed = _Guarded(workdir=str(tmpdir), path_state_file=state_file)
    resumed.system = _FakeSystem(par_file="p.yaml", outcomes=[(True, "")])
    resumed.run()

    assert(seen == [False, True])


def _fake_modules(system):
    class _Module:
        def check(self):
            pass
    return Dict(system=system, solver=_Module(), preprocess=_Module(),
                optimize=_Module())


def test_check_accepts_nonblocking_for_supported_workflow(tmpdir):
    workflow = Forward(workdir=str(tmpdir), modules=_fake_modules(
        _FakeSystem(par_file="p.yaml", blocking=False)))
    workflow.check()


def test_check_refuses_nonblocking_for_unsupported_workflow(tmpdir):
    """
    `Inversion`'s line search still blocks (and is recursive), so it must
    refuse a non-blocking system instead of half-working under it
    """
    workflow = Inversion(workdir=str(tmpdir), modules=_fake_modules(
        _FakeSystem(par_file="p.yaml", blocking=False)))
    with pytest.raises(AssertionError, match="does not yet support"):
        workflow.check()
