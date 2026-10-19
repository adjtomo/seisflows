"""
Test the TestFlow workflow, which is for finding out if a system module works
on a cluster. Here the cluster is simulated: a stand-in system runs the tasks
that are given to it (in this process, where a real one would on other nodes),
and keeps track of which failed for the workflow to ask about. No scheduler and
no SPECFEM is needed.
"""
import io
import os
import sys
import time
import pytest
import contextlib

from seisflows.tools.config import Dict, set_task_id
from seisflows.workflow.dispatch import WorkflowSuspended
from seisflows.workflow.test_flow import TestFlow


class _SimulatedScheduler:
    """
    A system module that runs tasks as soon as they are submitted, writing the
    output of each task to its own log file like a scheduler does. Mimics the
    blocking and non-blocking behavior of `system.Slurm`.
    """
    def __init__(self, workdir, blocking):
        self.blocking = blocking
        self.path = Dict(log_files=os.path.join(workdir, "logs"),
                         par_file=os.path.join(workdir, "parameters.yaml"))
        self.ntask = self.nproc = self.tasktime = self.walltime = None
        self.rerun = 0
        self.job_failures = {}  # job id: task ids that failed
        self.jobs = []  # (job id, task ids) of every job run
        self.submits = []  # what each resume job waits for
        self.times = (4., 62.)  # what the queue says

    def check(self):
        pass

    def setup(self):
        os.makedirs(self.path.log_files, exist_ok=True)

    def run(self, funcs, single=False, array=None, blocking=True, _attempts=0,
            **kwargs):
        if single:
            task_ids = [0]
        elif array:
            task_ids = [int(i) for i in str(array).split(",")]
        else:
            task_ids = list(range(self.ntask))

        job_id = str(100 + len(self.jobs))
        self.jobs.append((job_id, task_ids))
        failed = []
        for task_id in task_ids:
            set_task_id(task_id)
            output = io.StringIO()
            try:
                with contextlib.redirect_stdout(output):
                    for func in funcs:
                        func(**kwargs)
            except SystemExit:
                failed.append(task_id)
            with open(os.path.join(self.path.log_files,
                                   f"{job_id}_{task_id}"), "w") as f:
                f.write(output.getvalue())
        self.job_failures[job_id] = failed

        if not blocking:
            return job_id
        # The way that a blocking system recovers from failures
        if failed:
            if _attempts < self.rerun:
                return self.run(funcs, array=",".join(map(str, failed)),
                                _attempts=_attempts + 1, **kwargs)
            sys.exit(-1)

    def submit(self, workdir=None, parameter_file=None, depends_on=None):
        self.submits.append(depends_on)
        return f"resume-{depends_on}"

    def job_outcome(self, job_id):
        failed = self.job_failures[job_id]
        return not failed, ",".join(str(i) for i in failed)

    def job_times(self, job_id):
        return self.times


@pytest.fixture(autouse=True)
def _no_waiting(monkeypatch):
    """The tests wait for queues to catch up, which a simulation doesn't need"""
    monkeypatch.setattr(time, "sleep", lambda seconds: None)


def _testflow(tmpdir, system):
    """A test workflow as one process of a workflow would make it"""
    workflow = TestFlow(modules=Dict(system=system, solver=None,
                                     preprocess=None, optimize=None),
                        workdir=str(tmpdir))
    workflow.check()
    workflow.setup()
    return workflow


def _run_to_completion(tmpdir, system, max_processes=20):
    """
    Run the workflow, starting a new one each time that one suspends to wait
    for a job like the scheduler would

    :rtype: (int, str)
    :return: the number of processes it took, and the final status
    """
    for process in range(1, max_processes + 1):
        status = _testflow(tmpdir, system).run()
        if status != "SUSPENDED":
            return process, status
    raise AssertionError(f"still not done after {max_processes} processes")


def test_blocking_testflow_passes_in_one_process(tmpdir):
    """The tests as they have always been, on a system that waits for jobs"""
    system = _SimulatedScheduler(str(tmpdir), blocking=True)
    processes, status = _run_to_completion(tmpdir, system)

    assert((processes, status) == (1, "COMPLETE"))
    assert(system.submits == [])
    # 3 and 1 and 3 tasks, then 3 and the 2 that failed, run again
    assert([len(ids) for _, ids in system.jobs] == [3, 1, 3, 3, 2])


def test_nonblocking_testflow_resumes_after_every_job(tmpdir):
    """
    Each test dispatches a job and the main job stops, to be resumed by the
    scheduler after it. Six jobs (the failed tasks of the re-run test are one
    job of their own) so six resumes, and the workflow finishes in the seventh.
    """
    system = _SimulatedScheduler(str(tmpdir), blocking=False)
    processes, status = _run_to_completion(tmpdir, system)

    assert((processes, status) == (7, "COMPLETE"))
    assert([len(ids) for _, ids in system.jobs] == [3, 1, 3, 3, 2, 1])
    # The main job waits on each job in turn
    assert(system.submits == [job_id for job_id, _ in system.jobs])


def test_nonblocking_failures_are_found_by_the_main_job(tmpdir):
    """The tests whose jobs fail on purpose were caught, not just missed"""
    system = _SimulatedScheduler(str(tmpdir), blocking=False)
    _run_to_completion(tmpdir, system)

    # single job, tasks 1 and 2 of three, then those two again
    assert([system.job_failures[job_id] for job_id, _ in system.jobs] ==
           [[], [0], [1, 2], [1, 2], [1, 2], [0]])


def test_nonblocking_testflow_runs_again(tmpdir):
    """After passing it starts over instead of skipping every finished test"""
    system = _SimulatedScheduler(str(tmpdir), blocking=False)
    _run_to_completion(tmpdir, system)
    njobs = len(system.jobs)

    processes, status = _run_to_completion(tmpdir, system)
    assert((processes, status) == (7, "COMPLETE"))
    assert(len(system.jobs) == 2 * njobs)


def test_preflight_test_is_only_for_nonblocking_jobs(tmpdir):
    blocking = _testflow(tmpdir, _SimulatedScheduler(str(tmpdir), True))
    nonblocking = _testflow(tmpdir, _SimulatedScheduler(str(tmpdir), False))
    names = [task.__name__ for task in nonblocking.task_list]

    assert("test_preflight_failure_is_not_retried" in names)
    assert("test_preflight_failure_is_not_retried" not in
           [task.__name__ for task in blocking.task_list])


def test_a_test_that_does_not_fail_when_it_should_is_an_error(tmpdir):
    """Used to be a log message that the workflow carried on after"""
    system = _SimulatedScheduler(str(tmpdir), blocking=False)
    system.job_outcome = lambda job_id: (True, "")  # lies that nothing failed

    with pytest.raises(AssertionError, match="expected to fail"):
        _run_to_completion(tmpdir, system)


def test_failed_check_on_a_node_is_not_retried_even_with_reruns(tmpdir):
    """Only one job runs, though two re-runs are allowed"""
    system = _SimulatedScheduler(str(tmpdir), blocking=False)
    _run_to_completion(tmpdir, system)

    preflight_jobs = [ids for _, ids in system.jobs][5:]
    assert(preflight_jobs == [[0]])
    assert(system.rerun == 2)


def test_missing_job_times_is_a_warning_not_a_failure(tmpdir, caplog):
    """The times are only for the log, a system that can't say still passes"""
    system = _SimulatedScheduler(str(tmpdir), blocking=False)
    system.times = None
    _, status = _run_to_completion(tmpdir, system)

    assert(status == "COMPLETE")
    assert("`system.job_times()` does not work" in caplog.text)


def test_job_times_are_reported(tmpdir, caplog):
    import logging
    caplog.set_level(logging.INFO, logger="seisflows")
    system = _SimulatedScheduler(str(tmpdir), blocking=False)
    _run_to_completion(tmpdir, system)

    assert("queued for 4s and ran for 1m02s" in caplog.text)


def test_testflow_is_not_collected_as_a_test():
    assert(TestFlow.__test__ is False)


def test_nothing_but_parameters_is_public(tmpdir):
    """
    `seisflows configure` writes every public attribute of a workflow to the
    parameter file, so objects that are not parameters must be private
    """
    workflow = TestFlow(workdir=str(tmpdir))
    assert(sorted(key for key in vars(workflow) if not key.startswith("_")) ==
           sorted(["stop_after", "generate_data", "export_traces",
                   "export_residuals", "path"]))
