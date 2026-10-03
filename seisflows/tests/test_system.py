"""
Test the system modules ability to run jobs
"""
import os
import pytest
import subprocess
from glob import glob
from unittest.mock import patch, MagicMock
from seisflows.tools import unix
from seisflows.tools.config import get_task_id
from seisflows.system.workstation import Workstation
from seisflows.system.cluster import Cluster
from seisflows.system.chinook import Chinook
from seisflows.system.slurm import Slurm
from seisflows.system.pbs import Pbs
from seisflows.system.miyabi import Miyabi


def _test_function_a(tmpdir, **kwargs):
    """
    A test function that simply prints, used to check run

    .. note::
        Must accept `**kwargs` because `Cluster.run()` pickles `verbose` and
        `log_level` alongside any user-provided kwargs and passes the full
        set to every function in the run list (see `system.cluster.Cluster.
        run()` and `system/runscripts/run`)
    """
    test_dir = os.path.join(tmpdir, "test_dir")
    if not os.path.exists(test_dir):
        unix.mkdir(test_dir)


def _test_function_b(tmpdir, **kwargs):
    """
    A test function that simply prints
    """
    idx = get_task_id()
    test_fid = f"test_file_{idx:0>2}.txt"
    test_file = os.path.join(tmpdir, "test_dir", test_fid)
    print(f"making test file {test_file}")
    open(test_file, "w").close()


def test_workstation_run(tmpdir, ntask=26):
    """
    Make sure that Workstation run simply calls the functions in question
    """
    system = Workstation(workdir=tmpdir, path_system=tmpdir, ntask=ntask)
    system.setup()  # make required dir. structure
    system.run(funcs=[_test_function_a, _test_function_b], tmpdir=tmpdir)

    assert(len(glob(os.path.join(system.path.log_files, "*"))) == ntask)
    assert(len(glob(os.path.join(tmpdir, "test_dir", "*"))) == ntask)


def test_workstation_run_single(tmpdir, ntask=26):
    """
    Make sure that Workstation run simply calls the functions only once
    """
    system = Workstation(workdir=tmpdir, path_system=tmpdir, ntask=ntask)
    system.setup()  # make required dir. structure
    system.run(funcs=[_test_function_a, _test_function_b], tmpdir=tmpdir,
               single=True)

    assert(len(glob(os.path.join(system.path.log_files, "*"))) == 1)
    assert(len(glob(os.path.join(tmpdir, "test_dir", "*"))) == 1)


def test_cluster_run(tmpdir, ntask=26):
    """
    Make sure that Cluster run can pickle the tasks and run them with subprocess
    """
    system = Cluster(workdir=tmpdir, path_system=tmpdir, ntask=ntask)
    system.setup()  # make required dir. structure
    system.run(funcs=[_test_function_a, _test_function_b], tmpdir=tmpdir)

    assert(len(glob(os.path.join(system.path.log_files, "*"))) == ntask)
    assert(len(glob(os.path.join(tmpdir, "test_dir", "*"))) == ntask)


def test_cluster_run_single(tmpdir, ntask=26):
    """
    Make sure that Cluster runs single tasks accetpably
    """
    system = Cluster(workdir=tmpdir, path_system=tmpdir, ntask=ntask)
    system.setup()  # make required dir. structure
    system.run(funcs=[_test_function_a, _test_function_b], tmpdir=tmpdir,
               single=True)

    assert(len(glob(os.path.join(system.path.log_files, "*"))) == 1)
    assert(len(glob(os.path.join(tmpdir, "test_dir", "*"))) == 1)


# -----------------------------------------------------------------------------
# Non-blocking master job: scheduler-native dependency-chaining primitives.
# These mock `subprocess.run` rather than touch a real scheduler -- see
# `seisflows/tests/conftest.py` for why that matters here in particular
# (`Cluster.run()`/`Slurm.run()` spawn real subprocesses for their *compute*
# tasks, but the dependency-chaining calls tested below go through
# `subprocess.run` directly in-process, so they don't hit that gotcha; mocking
# keeps these tests fast and independent of any real cluster regardless).
# -----------------------------------------------------------------------------

def test_slurm_dependency_flag(tmpdir):
    """
    `Slurm._dependency_flag()` should build SLURM's 'afterany' dependency
    syntax, which releases the dependent submission once the given job ID
    finishes, successfully or not
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1)
    assert(system._dependency_flag("123456") == "--dependency=afterany:123456")


def test_cluster_dependency_flag_not_implemented(tmpdir):
    """
    Generalized `Cluster` mimics a scheduler locally (see module docstring)
    and has no real external job to depend on, so `_dependency_flag()` (and
    by extension `submit(depends_on=...)`) should refuse rather than build
    a meaningless call
    """
    system = Cluster(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1)
    with pytest.raises(NotImplementedError):
        system._dependency_flag("123456")


def test_submit_depends_on_appends_dependency_flag(tmpdir):
    """
    `submit(depends_on=job_id)` should append the system's dependency flag
    to the submit call's header, so the scheduler won't release the
    resubmitted master job until `job_id` finishes
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1,
                     submit_to="debug")
    system.setup()

    with patch("seisflows.system.cluster.subprocess.run") as mock_run:
        system.submit(workdir=str(tmpdir), parameter_file="parameters.yaml",
                      depends_on="123456")

    assert(mock_run.called)
    submit_call = mock_run.call_args[0][0]
    assert("--dependency=afterany:123456" in submit_call)
    assert("--partition=debug" in submit_call)  # submit_to, not partition
    assert(system.submit_workflow in submit_call)


def test_submit_without_depends_on_omits_dependency_flag(tmpdir):
    """A normal (first) submission should not carry a dependency flag"""
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1,
                     submit_to="debug")
    system.setup()

    with patch("seisflows.system.cluster.subprocess.run") as mock_run:
        system.submit(workdir=str(tmpdir), parameter_file="parameters.yaml")

    submit_call = mock_run.call_args[0][0]
    assert("--dependency" not in submit_call)


def test_submit_direct_ignores_depends_on(tmpdir):
    """
    `direct=True` runs the master job immediately in-process (see
    `Cluster.submit()`), so there is no scheduler call to attach a
    dependency to -- `depends_on` should be silently ignored, not raise
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1,
                     submit_to="debug")
    system.setup()

    with patch("seisflows.system.cluster.subprocess.run") as mock_run:
        system.submit(workdir=str(tmpdir), parameter_file="parameters.yaml",
                      direct=True, depends_on="123456")

    submit_call = mock_run.call_args[0][0]
    assert("--dependency" not in submit_call)


def test_slurm_run_blocking_false_returns_job_id_without_monitoring(tmpdir):
    """
    `run(blocking=False)` should submit the job and return its ID
    immediately, without calling `monitor_job_status` (i.e., without
    polling/waiting for it to finish)
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=2,
                     submit_to="debug")
    system.setup()

    mock_submit_result = MagicMock()
    mock_submit_result.stdout = "654321\n"

    with patch("seisflows.system.slurm.subprocess.run",
              return_value=mock_submit_result) as mock_run, \
         patch.object(system, "monitor_job_status") as mock_monitor:
        job_id = system.run(funcs=[_test_function_a], tmpdir=str(tmpdir),
                            blocking=False)

    assert(job_id == "654321")
    mock_monitor.assert_not_called()
    assert(mock_run.called)


def test_slurm_run_blocking_true_still_monitors(tmpdir):
    """
    Default (`blocking=True`) behavior must be unchanged: `run()` still
    calls `monitor_job_status` and waits, exactly as before this feature
    existed
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=2,
                     submit_to="debug")
    system.setup()

    mock_submit_result = MagicMock()
    mock_submit_result.stdout = "654321\n"

    with patch("seisflows.system.slurm.subprocess.run",
              return_value=mock_submit_result), \
         patch.object(system, "monitor_job_status",
                      return_value=1) as mock_monitor:
        job_id = system.run(funcs=[_test_function_a], tmpdir=str(tmpdir))

    mock_monitor.assert_called_once_with("654321")
    assert(job_id is None)  # blocking run has no meaningful return value




def test_submit_resume_appends_log_and_skips_backup(tmpdir):
    """
    A dependent (resumed) submission should make the scheduler APPEND to the
    output log rather than truncate it, and so should NOT back up the log
    first (backing up on every resume would pile up numbered copies)
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1,
                     submit_to="debug")
    system.setup()
    with open(system.path.output_log, "w") as f:
        f.write("log from an earlier master job invocation\n")

    with patch("seisflows.system.cluster.subprocess.run") as mock_run:
        system.submit(workdir=str(tmpdir), parameter_file="parameters.yaml",
                      depends_on="123456")

    assert("--open-mode=append" in mock_run.call_args[0][0])
    assert(glob(os.path.join(system.path.log_files, "*")) == [])


def test_submit_first_submission_backs_up_log(tmpdir, monkeypatch):
    """A normal (first) submission still backs up an existing log"""
    # `copy_file()` resolves the file relative to the cwd, which is always
    # the workdir when SeisFlows runs for real
    monkeypatch.chdir(tmpdir)
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1,
                     submit_to="debug")
    system.setup()
    with open(system.path.output_log, "w") as f:
        f.write("log from a previous workflow\n")

    with patch("seisflows.system.cluster.subprocess.run") as mock_run:
        system.submit(workdir=str(tmpdir), parameter_file="parameters.yaml")

    assert("--open-mode" not in mock_run.call_args[0][0])
    assert(len(glob(os.path.join(system.path.log_files, "*"))) == 1)


def test_submit_resume_failure_is_not_silent(tmpdir):
    """
    If the dependent job fails to submit, the workflow would exit after
    dispatching and never be resumed -- that must be a hard error, not
    ignored like a failed first submission is
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1,
                     submit_to="debug")
    system.setup()

    error = subprocess.CalledProcessError(1, "sbatch")
    with patch("seisflows.system.cluster.subprocess.run", side_effect=error):
        with pytest.raises(SystemExit):
            system.submit(workdir=str(tmpdir), depends_on="123456")


def test_slurm_job_outcome_success(tmpdir):
    """All tasks completed: success, nothing to re-run"""
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=2)
    with patch.object(system, "monitor_job_status", return_value=1):
        assert(system.job_outcome("99") == (True, ""))


def test_slurm_job_outcome_reports_failed_task_ids(tmpdir):
    """
    Failed task ids come from each array element's own id ('99_3' -> '3')
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=4)
    jobs = ("99_0", "99_1", "99_2", "99_3")
    states = ("COMPLETED", "FAILED", "COMPLETED", "TIMEOUT")
    with patch.object(system, "monitor_job_status", return_value=-1), \
         patch.object(system, "query_job_states",
                      return_value=(jobs, states)):
        assert(system.job_outcome("99") == (False, "1,3"))


def test_slurm_job_outcome_after_partial_rerun(tmpdir):
    """
    A re-run of only tasks 1 and 3 produces an array with just those two
    elements. If task 3 fails again, it must be reported as '3' -- reading
    the position in the list instead would wrongly give '1'
    """
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=4)
    jobs = ("100_1", "100_3")
    states = ("COMPLETED", "FAILED")
    with patch.object(system, "monitor_job_status", return_value=-1), \
         patch.object(system, "query_job_states",
                      return_value=(jobs, states)):
        assert(system.job_outcome("100") == (False, "3"))


def test_cluster_job_outcome_not_implemented(tmpdir):
    """Generalized `Cluster` has no real scheduler to inspect"""
    system = Cluster(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1)
    with pytest.raises(NotImplementedError):
        system.job_outcome("99")


def test_nonblocking_support_flags():
    """
    Support must be declared per class rather than inferred: only systems
    that implement `run(blocking=False)` and `job_outcome()` claim it
    """
    assert(Cluster._supports_nonblocking is False)
    assert(Slurm._supports_nonblocking is True)
    assert(Chinook._supports_nonblocking is True)


def test_check_refuses_nonblocking_on_unsupported_system(tmpdir):
    """`blocking=False` on a system that can't honor it fails at check()"""
    system = Cluster(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1,
                     blocking=False)
    with patch("seisflows.system.workstation.Workstation.check"):
        with pytest.raises(AssertionError, match="non-blocking"):
            system.check()


def test_check_allows_nonblocking_on_supported_system(tmpdir):
    """...and passes for one that can"""
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=1,
                     blocking=False)
    with patch("seisflows.system.workstation.Workstation.check"):
        system.check()


# -----------------------------------------------------------------------------
# Non-blocking master job on PBS (Miyabi). As above, `subprocess.run` is mocked
# -----------------------------------------------------------------------------

def _miyabi(tmpdir, ntask=4, **kwargs):
    """Miyabi system, set up, with the minimum parameters to submit"""
    system = Miyabi(workdir=str(tmpdir), path_system=str(tmpdir),
                    ntask=ntask, group="gz00", **kwargs)
    system.setup()
    return system


def _qsub_result(*job_ids):
    """Mock `subprocess.run` results, one per qsub, each a job id"""
    results = []
    for job_id in job_ids:
        result = MagicMock()
        result.stdout = f"{job_id}\n"
        results.append(result)
    return results


def test_pbs_dependency_flag(tmpdir):
    """
    PBS's 'afterany' dependency, with the ',' that `run(blocking=False)`
    joins several job ids with turned into PBS's ':'
    """
    system = _miyabi(tmpdir)
    assert(system._dependency_flag("123.opbs") ==
           "-W depend=afterany:123.opbs")
    assert(system._dependency_flag("123[].opbs,124[].opbs") ==
           "-W depend=afterany:123[].opbs:124[].opbs")


def test_pbs_submit_depends_on_appends_dependency_flag(tmpdir):
    """A resumed master job is submitted with the dependency, a first one
    without"""
    system = _miyabi(tmpdir)

    with patch("seisflows.system.pbs.subprocess.run") as mock_run:
        system.submit(workdir=str(tmpdir), depends_on="123.opbs")
    assert("-W depend=afterany:123.opbs" in mock_run.call_args[0][0])

    with patch("seisflows.system.pbs.subprocess.run") as mock_run:
        system.submit(workdir=str(tmpdir))
    assert("-W depend=" not in mock_run.call_args[0][0])


def test_pbs_submit_resume_failure_is_not_silent(tmpdir):
    """A failed submission of the resume job must stop the workflow with an
    error that says how to recover"""
    system = _miyabi(tmpdir)
    error = subprocess.CalledProcessError(1, "qsub")
    with patch("seisflows.system.pbs.subprocess.run", side_effect=error), \
         patch("seisflows.system.pbs.logger") as mock_logger:
        with pytest.raises(SystemExit):
            system.submit(workdir=str(tmpdir), depends_on="123.opbs")
    assert("resume manually" in mock_logger.critical.call_args[0][0])


def test_pbs_run_blocking_false_returns_job_id_without_monitoring(tmpdir):
    """Submits an array job and returns at once, with the id of the array"""
    system = _miyabi(tmpdir)
    with patch("seisflows.system.pbs.subprocess.run",
               side_effect=_qsub_result("55[].opbs")) as mock_run, \
         patch.object(system, "monitor_job_status") as mock_monitor:
        job_id = system.run(funcs=[_test_function_a], tmpdir=str(tmpdir),
                            blocking=False)

    assert(job_id == "55[].opbs")
    mock_monitor.assert_not_called()
    assert(mock_run.call_count == 1)
    assert("-J 0-3%" in mock_run.call_args[0][0])


def test_pbs_run_blocking_true_still_monitors(tmpdir):
    """The default path is unchanged by the refactor into
    `_submit_run_call`"""
    system = _miyabi(tmpdir)
    with patch("seisflows.system.pbs.subprocess.run",
               side_effect=_qsub_result("55[].opbs")), \
         patch("seisflows.system.pbs.time.sleep"), \
         patch.object(system, "monitor_job_status",
                      return_value=1) as mock_monitor:
        assert(system.run(funcs=[_test_function_a], tmpdir=str(tmpdir))
               is None)

    mock_monitor.assert_called_once_with("55[].opbs")


def test_pbs_run_rerun_submits_one_job_per_contiguous_range(tmpdir):
    """
    PBS's '-J' takes a single contiguous range, so re-running failed tasks
    1,3,4 (as the workflow asks for them) takes two jobs, 1-1 and 3-4, whose
    ids are returned together
    """
    system = _miyabi(tmpdir, ntask=6)
    with patch("seisflows.system.pbs.subprocess.run",
               side_effect=_qsub_result("56[].opbs", "57[].opbs")) as mock_run:
        job_id = system.run(funcs=[_test_function_a], tmpdir=str(tmpdir),
                            blocking=False, array="1,3,4")

    assert(job_id == "56[].opbs,57[].opbs")
    calls = [c[0][0] for c in mock_run.call_args_list]
    assert(len(calls) == 2)
    assert("-J 1-1 " in calls[0])
    assert("-J 3-4 " in calls[1])


def test_pbs_run_single_task_is_not_an_array(tmpdir):
    """A single-task job has no '-J', so a re-run of task 0 submits as is"""
    system = _miyabi(tmpdir, ntask=1)
    with patch("seisflows.system.pbs.subprocess.run",
               side_effect=_qsub_result("58.opbs")) as mock_run:
        job_id = system.run(funcs=[_test_function_a], tmpdir=str(tmpdir),
                            blocking=False, array="0")

    assert(job_id == "58.opbs")
    assert("-J" not in mock_run.call_args[0][0])


def test_pbs_array_indices():
    """Expands the index lists and ranges that failed tasks arrive as"""
    assert(Pbs._array_indices("1,3,5-7") == [1, 3, 5, 6, 7])
    assert(Pbs._array_indices("2") == [2])
    assert(Pbs._array_indices("4,2,2-3") == [2, 3, 4])


def test_pbs_job_outcome_success(tmpdir):
    """All tasks completed: success, nothing to re-run"""
    system = _miyabi(tmpdir)
    with patch.object(system, "monitor_job_status", return_value=1):
        assert(system.job_outcome("55[].opbs") == (True, ""))


def test_pbs_job_outcome_reports_failed_task_ids(tmpdir):
    """Failed task ids come from the subjob ids ('55[3].opbs' -> 3)"""
    system = _miyabi(tmpdir)
    jobs = ("55[0].opbs", "55[1].opbs", "55[2].opbs", "55[3].opbs")
    states = ("COMPLETED", "FAILED", "COMPLETED", "FAILED")
    with patch.object(system, "monitor_job_status", return_value=-1), \
         patch.object(system, "query_job_states",
                      return_value=(jobs, states)):
        assert(system.job_outcome("55[].opbs") == (False, "1,3"))


def test_pbs_job_outcome_across_several_jobs(tmpdir):
    """
    A re-run that took two jobs is judged as a whole: monitored together,
    and the failures of both are reported
    """
    system = _miyabi(tmpdir)
    answers = {"56[].opbs": (("56[1].opbs",), ("COMPLETED",)),
               "57[].opbs": (("57[3].opbs", "57[4].opbs"),
                             ("FAILED", "COMPLETED"))}
    with patch.object(system, "monitor_job_status",
                      return_value=-1) as mock_monitor, \
         patch.object(system, "query_job_states",
                      side_effect=lambda jid, sort=False: answers[jid]):
        outcome = system.job_outcome("56[].opbs,57[].opbs")

    assert(outcome == (False, "3"))
    assert(mock_monitor.call_args[0][0] == ["56[].opbs", "57[].opbs"])


def test_pbs_job_outcome_single_job_failure_is_task_zero(tmpdir):
    """A job that isn't an array has no index and is task 0"""
    system = _miyabi(tmpdir, ntask=1)
    with patch.object(system, "monitor_job_status", return_value=-1), \
         patch.object(system, "query_job_states",
                      return_value=(("58.opbs",), ("FAILED",))):
        assert(system.job_outcome("58.opbs") == (False, "0"))


def test_monitor_job_status_waits_for_every_job_in_a_list(tmpdir):
    """
    If the scheduler can't report on one of several jobs yet, the others
    finishing is not enough: that is not complete, and used to be a crash
    """
    system = _miyabi(tmpdir)
    answers = iter([(("1.opbs",), ("COMPLETED",)), (None, None),
                    (("1.opbs",), ("COMPLETED",)),
                    (("2.opbs",), ("COMPLETED",))])
    with patch("seisflows.system.cluster.time.sleep"), \
         patch.object(system, "query_job_states",
                      side_effect=lambda jid, sort=False: next(answers)):
        assert(system.monitor_job_status(["1.opbs", "2.opbs"]) == 1)


def test_pbs_nonblocking_support_flags(tmpdir):
    """PBS and its children claim support and pass `check()`"""
    assert(Pbs._supports_nonblocking is True)
    assert(Miyabi._supports_nonblocking is True)
    system = _miyabi(tmpdir, blocking=False)
    with patch("seisflows.system.workstation.Workstation.check"):
        system.check()
