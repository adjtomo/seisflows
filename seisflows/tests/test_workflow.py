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

from seisflows.workflow.forward import Forward
from seisflows.workflow.inversion import Inversion


class _DummyForward(Forward):
    """
    Minimal `Forward` subclass with a 3-task list of no-op functions, so that
    `Forward.run()`'s checkpoint/resume/failure logic can be tested without
    any real system/solver modules.
    """
    def __init__(self, fail_on=None, **kwargs):
        self.calls = []
        self._fail_on = fail_on or set()
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

    def task_c(self):
        self.calls.append("task_c")


def test_forward_run_executes_all_tasks_in_order(tmpdir):
    """A fresh workflow should run every task, in order, exactly once"""
    workflow = _DummyForward(workdir=str(tmpdir))
    workflow.path.state_file = os.path.join(str(tmpdir), "sfstate.txt")
    workflow.setup = lambda: None  # setup() needs real modules, skip it

    # `setup()` would normally create the state file; do the minimal
    # equivalent by hand
    with open(workflow.path.state_file, "w") as f:
        f.write("# SeisFlows State File\n")

    workflow.run()

    assert(workflow.calls == ["task_a", "task_b", "task_c"])
    assert(workflow._states == {"task_a": 1, "task_b": 1, "task_c": 1})


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
    workflow.run()

    assert(workflow.calls == ["task_a", "task_b"])
    assert(workflow._states["task_a"] == 1)
    assert(workflow._states["task_b"] == 1)
    assert(workflow._states["task_c"] == 0)


class _DummyInversion(Inversion):
    """
    Minimal `Inversion` subclass with a single dummy task, used to test the
    per-iteration bookkeeping in `Inversion.run()`/`Inversion.checkpoint()`
    without needing real system/solver/preprocess/optimize modules
    """
    def __init__(self, **kwargs):
        self.calls = []
        super().__init__(**kwargs)

    @property
    def task_list(self):
        return [self.only_task]

    def only_task(self):
        self.calls.append(self.iteration)


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
