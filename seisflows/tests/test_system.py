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
