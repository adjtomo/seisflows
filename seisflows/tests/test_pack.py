"""
Test pack mode: running the tasks of an array in fewer jobs
"""
import json
import os
import pytest
import subprocess
import sys
from unittest.mock import patch, MagicMock

from seisflows.system import pack
from seisflows.system.chinook import Chinook
from seisflows.system.pack import (parse_array, split_tasks, pack_tasktime,
                                   job_key, failed_tasks, run_chunk)
from seisflows.tools.config import pickle_function_list
from seisflows.tests.test_system import (_test_function_a, _test_function_b,
                                         _miyabi, _qsub_result)


def test_parse_array():
    assert(parse_array("0,2,5-7") == [0, 2, 5, 6, 7])
    assert(parse_array("0-79%15") == list(range(80)))  # throttle is ignored
    assert(parse_array("3") == [3])
    assert(parse_array("") == [])


@pytest.mark.parametrize("ntask,njob,sizes", [
    (54, 4, [14, 14, 13, 13]),
    (54, 3, [18, 18, 18]),
    (54, 1, [54]),
    (10, 4, [3, 3, 2, 2]),
    (5, 9, [1, 1, 1, 1, 1]),  # never more jobs than tasks
])
def test_split_tasks_is_even_and_contiguous(ntask, njob, sizes):
    blocks = split_tasks(list(range(ntask)), njob)
    assert([len(block) for block in blocks] == sizes)
    # Neighbouring ids stay together, every task once and in order
    assert(sum(blocks, []) == list(range(ntask)))


def test_split_tasks_of_a_rerun():
    """The ids that a rerun is given are not a range"""
    assert(split_tasks([1, 4, 5, 9, 12], 2) == [[1, 4, 5], [9, 12]])


def test_pack_tasktime_covers_every_task_and_a_buffer():
    assert(pack_tasktime(10, 18) == 198)  # 3 hours 18 min for 18 x 10 min
    assert(pack_tasktime(10, 1) == 11)
    assert(isinstance(pack_tasktime(0.5, 3), int))


def test_job_key_is_the_same_for_the_array_and_its_elements():
    for job_id in ["123", "123[].opbs", "123[4].opbs", "123_4"]:
        assert(job_key(job_id) == "123")


# -----------------------------------------------------------------------------
# Running the tasks of a block
# -----------------------------------------------------------------------------

pack_script = Chinook.run_chunk

STUB = '''
import os, sys, time
task = int(os.environ["SEISFLOWS_TASKID"])
print("task", task, "of", os.environ.get("SLURM_ARRAY_TASK_ID"))
if task == 4:
    sys.exit(3)
if task == 5:
    time.sleep(30)
'''


def _stub_run_script(tmpdir):
    """A stand-in for the `run` script, so that the tasks are quick"""
    path = os.path.join(str(tmpdir), "stub_run.py")
    with open(path, "w") as f:
        f.write(STUB)
    return path


def test_run_chunk_runs_every_task_and_marks_each(tmpdir, monkeypatch):
    """
    A task that fails, and one that runs out of time, do not stop the tasks
    after them. Each has a log and a marker. Every task gets its own id even
    though the scheduler gave the whole block one array index
    """
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "0")
    markers = os.path.join(str(tmpdir), "markers")
    logdir = os.path.join(str(tmpdir), "logs")

    exit_code = run_chunk([3, 4, 5, 6], funcs="f.p", kwargs="k.p",
                          tasktime=0.03, logdir=logdir, markers=markers,
                          job_id="77", run_script=_stub_run_script(tmpdir))

    assert(exit_code == 1)
    found = sorted(os.listdir(os.path.join(markers, "77")))
    assert(found == ["3.ok", "4.fail", "5.fail", "6.ok"])

    log = lambda task: open(os.path.join(logdir, f"77_{task}.log")).read()
    assert("task 3 of 0" in log(3))  # not the array index, the task id
    assert("exited with code 3" in log(4))
    assert("exceeded its time limit" in log(5))
    assert("task 6" in log(6))


def test_run_chunk_succeeds_when_every_task_does(tmpdir):
    markers = os.path.join(str(tmpdir), "markers")
    exit_code = run_chunk([0, 1, 2], funcs="f.p", kwargs="k.p", tasktime=1,
                          logdir=os.path.join(str(tmpdir), "logs"),
                          markers=markers, job_id="8[1].opbs",
                          run_script=_stub_run_script(tmpdir))
    assert(exit_code == 0)
    assert(sorted(os.listdir(os.path.join(markers, "8"))) ==
           ["0.ok", "1.ok", "2.ok"])


def test_run_chunk_with_the_real_run_script(tmpdir):
    """The tasks are run by the usual `run` script, by their own ids"""
    test_dir = os.path.join(str(tmpdir), "test_dir")
    os.makedirs(test_dir)
    funcs, kwargs = pickle_function_list(
        [_test_function_b], path=str(tmpdir), tmpdir=str(tmpdir),
        verbose=False, level="INFO")

    exit_code = run_chunk([2, 3, 7], funcs=funcs, kwargs=kwargs, tasktime=2,
                          logdir=os.path.join(str(tmpdir), "logs"),
                          markers=os.path.join(str(tmpdir), "markers"),
                          job_id="5")

    assert(exit_code == 0)
    assert(sorted(os.listdir(test_dir)) == ["test_file_02.txt",
                                            "test_file_03.txt",
                                            "test_file_07.txt"])


def test_run_chunk_script_runs_the_block_of_its_array_index(tmpdir):
    """
    What a job does: the array element picks its block from the file of
    blocks, and names its markers after the array job
    """
    test_dir = os.path.join(str(tmpdir), "test_dir")
    os.makedirs(test_dir)
    funcs, kwargs = pickle_function_list(
        [_test_function_b], path=str(tmpdir), tmpdir=str(tmpdir),
        verbose=False, level="INFO")
    chunks = os.path.join(str(tmpdir), "chunks.json")
    with open(chunks, "w") as f:
        json.dump({"0": [0, 1], "1": [2, 3, 4]}, f)

    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(pack.__file__))))
    env = dict(os.environ, SLURM_ARRAY_TASK_ID="1", SLURM_ARRAY_JOB_ID="55",
               PYTHONPATH=os.pathsep.join([root] + sys.path))
    # Another test may have set it, and it is read before the array index
    env.pop("SEISFLOWS_TASKID", None)
    result = subprocess.run(
        [sys.executable, pack_script, "--funcs", funcs, "--kwargs", kwargs,
         "--chunks", chunks, "--tasktime", "2", "--logdir",
         os.path.join(str(tmpdir), "logs"), "--markers",
         os.path.join(str(tmpdir), "markers"), "--environment"],
        env=env, capture_output=True, text=True)

    assert(result.returncode == 0), result.stderr
    assert(sorted(os.listdir(test_dir)) == ["test_file_02.txt",
                                            "test_file_03.txt",
                                            "test_file_04.txt"])
    assert(sorted(os.listdir(os.path.join(str(tmpdir), "markers", "55"))) ==
           ["2.ok", "3.ok", "4.ok"])


def test_failed_tasks(tmpdir):
    markers = str(tmpdir)
    assert(failed_tasks(markers, "9") is None)  # not a packed job

    path = os.path.join(markers, "9")
    os.makedirs(path)
    with open(os.path.join(path, "chunks.json"), "w") as f:
        json.dump({"0": [0, 1, 2], "1": [3, 4]}, f)
    for marker in ["0.ok", "1.fail", "2.ok", "3.ok"]:
        open(os.path.join(path, marker), "w").close()

    # 1 failed and 4 never ran, e.g., the job hit its time limit
    assert(failed_tasks(markers, "9[].opbs") == [1, 4])


# -----------------------------------------------------------------------------
# Slurm submits packed jobs. `sbatch` is mocked
# -----------------------------------------------------------------------------

def _packed_chinook(tmpdir, ntask=6, njob=2, **kwargs):
    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir),
                     ntask=ntask, njob=njob, pack=True, tasktime=10, nproc=4,
                     **kwargs)
    system.setup()
    return system


def _sbatch(*job_ids):
    results = []
    for job_id in job_ids:
        result = MagicMock()
        result.stdout = f"{job_id}\n"
        results.append(result)
    return results


def test_slurm_run_packs_the_tasks_into_njob_jobs(tmpdir):
    system = _packed_chinook(tmpdir, ntask=7, njob=2)

    with patch("seisflows.system.slurm.subprocess.run",
               side_effect=_sbatch("321")) as mock_run:
        job_id = system.run([_test_function_a], tmpdir=str(tmpdir),
                            blocking=False)

    assert(job_id == "321")
    assert(mock_run.call_count == 1)
    call = mock_run.call_args[0][0]
    assert("--array=0-1 " in call)  # two jobs, not seven
    assert(system.run_chunk in call)
    assert("--tasktime 10 " in call)  # a task's limit, not the job's
    assert("--time=44 " in call)  # 4 tasks x 10 min + 10%

    # The blocks that the jobs run, for the jobs and for `job_outcome`
    chunks_fid = [word for word in call.split() if word.endswith("chunks.json")]
    for fid in [chunks_fid[0],
                os.path.join(system._markers, "321", "chunks.json")]:
        with open(fid) as f:
            assert(json.load(f) == {"0": [0, 1, 2, 3], "1": [4, 5, 6]})


def test_slurm_run_is_the_usual_array_when_not_packed(tmpdir):
    system = _packed_chinook(tmpdir, ntask=4, njob=2)
    system.pack = False
    with patch("seisflows.system.slurm.subprocess.run",
               side_effect=_sbatch("321")) as mock_run:
        system.run([_test_function_a], tmpdir=str(tmpdir), blocking=False)
    call = mock_run.call_args[0][0]
    assert("--array=0-3%" in call)
    assert(system.run_chunk not in call)
    assert("--time=10 " in call)


def test_slurm_run_does_not_pack_a_rerun_of_a_few_tasks(tmpdir):
    """Two failed tasks and `njob` of 2: a job each, the usual array"""
    system = _packed_chinook(tmpdir, ntask=54, njob=2)
    with patch("seisflows.system.slurm.subprocess.run",
               side_effect=_sbatch("322")) as mock_run:
        system.run([_test_function_a], tmpdir=str(tmpdir), blocking=False,
                   array="4,9")
    call = mock_run.call_args[0][0]
    assert("--array=4,9 " in call)
    assert(system.run_chunk not in call)
    assert("--time=10 " in call)


def test_slurm_run_packs_a_rerun_of_many_tasks(tmpdir):
    system = _packed_chinook(tmpdir, ntask=54, njob=2)
    with patch("seisflows.system.slurm.subprocess.run",
               side_effect=_sbatch("322")) as mock_run:
        system.run([_test_function_a], tmpdir=str(tmpdir), blocking=False,
                   array="1,4-6,9")
    call = mock_run.call_args[0][0]
    assert("--array=0-1 " in call)
    with open(os.path.join(system._markers, "322", "chunks.json")) as f:
        assert(json.load(f) == {"0": [1, 4, 5], "1": [6, 9]})


def test_slurm_single_job_is_never_packed(tmpdir):
    system = _packed_chinook(tmpdir, ntask=6, njob=2)
    with patch("seisflows.system.slurm.subprocess.run",
               side_effect=_sbatch("323")) as mock_run:
        system.run([_test_function_a], tmpdir=str(tmpdir), single=True,
                   blocking=False)
    assert(system.run_chunk not in mock_run.call_args[0][0])


def test_slurm_job_outcome_of_a_packed_job(tmpdir):
    """It is the tasks that failed, whatever state the scheduler says"""
    system = _packed_chinook(tmpdir)
    path = os.path.join(system._markers, "99")
    os.makedirs(path)
    with open(os.path.join(path, "chunks.json"), "w") as f:
        json.dump({"0": [0, 1, 2], "1": [3, 4, 5]}, f)
    for marker in ["0.ok", "1.ok", "2.ok", "3.ok", "4.fail"]:  # 5 never ran
        open(os.path.join(path, marker), "w").close()

    with patch.object(system, "monitor_job_status") as mock_monitor:
        assert(system.job_outcome("99") == (False, "4,5"))
    mock_monitor.assert_not_called()

    open(os.path.join(path, "4.ok"), "w").close()
    open(os.path.join(path, "5.ok"), "w").close()
    assert(system.job_outcome("99") == (True, ""))


def test_slurm_blocking_run_reruns_the_tasks_that_failed(tmpdir):
    """The blocking path finds the failed tasks from the markers too"""
    system = _packed_chinook(tmpdir, ntask=6, njob=3, rerun=1)
    calls = []

    def sbatch(call, **kwargs):
        calls.append(call)
        if len(calls) == 1:
            # What the three jobs do: task 4 fails, so does the job that ran it
            path = os.path.join(system._markers, "400")
            os.makedirs(path)
            for marker in ["0.ok", "1.ok", "2.ok", "3.ok", "4.fail", "5.ok"]:
                open(os.path.join(path, marker), "w").close()
            return _sbatch("400")[0]
        return _sbatch("401")[0]

    with patch("seisflows.system.slurm.subprocess.run", side_effect=sbatch), \
         patch("seisflows.system.slurm.time.sleep"), \
         patch.object(system, "monitor_job_status", side_effect=[-1, 1]):
        system.run([_test_function_a], tmpdir=str(tmpdir))

    assert(len(calls) == 2)
    assert("--array=4 " in calls[1])
    assert(system.run_chunk not in calls[1])


def test_check_pack_needs_a_valid_njob(tmpdir):
    with open(os.path.join(str(tmpdir), "parameters.yaml"), "w") as f:
        f.write("# parameters\n")
    for njob in [None, 0, 7, 1.5]:
        system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir),
                         ntask=6, pack=True, njob=njob)
        with pytest.raises(AssertionError) as e:
            system.check()
        assert("njob" in str(e.value))

    system = Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=6,
                     pack=True, njob=3)
    system.check()
    # Ignored unless packing
    Chinook(workdir=str(tmpdir), path_system=str(tmpdir), ntask=6,
            njob=None).check()


# -----------------------------------------------------------------------------
# PBS submits packed jobs. `qsub` is mocked
# -----------------------------------------------------------------------------

def test_pbs_run_packs_the_tasks_into_njob_jobs(tmpdir):
    system = _miyabi(tmpdir, ntask=54, njob=3, pack=True, tasktime=10,
                     partition_compute="short-g")

    with patch("seisflows.system.pbs.subprocess.run",
               side_effect=_qsub_result("500[].opbs")) as mock_run:
        job_id = system.run([_test_function_a], tmpdir=str(tmpdir),
                            blocking=False)

    assert(job_id == "500[].opbs")
    assert(mock_run.call_count == 1)  # one array, not a job per range
    call = mock_run.call_args[0][0]
    assert(" -J 0-2 " in call)
    assert(" -q short-g " in call)
    assert(" -l walltime=03:18:00 " in call)  # 18 tasks x 10 min + 10%
    assert("SEISFLOWS_ENTRY_POINT=run_chunk," in call)
    assert("SEISFLOWS_TASKTIME=10," in call)
    assert(f"SEISFLOWS_MARKERS={system._markers}" in call)
    assert(" -v " in call and "SEISFLOWS_TASKID" not in call)

    with open(os.path.join(system._markers, "500", "chunks.json")) as f:
        blocks = json.load(f)
    assert([len(block) for block in blocks.values()] == [18, 18, 18])
    assert(blocks["1"][0] == 18)


def test_pbs_run_of_one_packed_job_is_not_an_array(tmpdir):
    """Some PBS implementations reject an array range of one, '0-0', so a
    single job runs every task as the block with index 0"""
    system = _miyabi(tmpdir, ntask=6, njob=1, pack=True, tasktime=10)

    with patch("seisflows.system.pbs.subprocess.run",
               side_effect=_qsub_result("501.opbs")) as mock_run:
        job_id = system.run([_test_function_a], tmpdir=str(tmpdir),
                            blocking=False)

    call = mock_run.call_args[0][0]
    assert(" -J " not in call)
    assert("SEISFLOWS_TASKID=0" in call)
    assert("SEISFLOWS_ENTRY_POINT=run_chunk," in call)
    assert(" -l walltime=01:06:00 " in call)  # 6 tasks x 10 min + 10%
    with open(os.path.join(system._markers, "501", "chunks.json")) as f:
        assert(json.load(f) == {"0": [0, 1, 2, 3, 4, 5]})
    assert(job_id == "501.opbs")


def test_pbs_run_is_the_usual_array_when_not_packed(tmpdir):
    system = _miyabi(tmpdir, ntask=6, njob=2, tasktime=10)  # `pack` is False
    with patch("seisflows.system.pbs.subprocess.run",
               side_effect=_qsub_result("502[].opbs")) as mock_run:
        system.run([_test_function_a], tmpdir=str(tmpdir), blocking=False)
    call = mock_run.call_args[0][0]
    assert("SEISFLOWS_ENTRY_POINT=run," in call)
    assert(" -J 0-5%" in call)


def test_pbs_job_outcome_of_a_packed_job(tmpdir):
    """The job id has changed form since it was submitted, the array
    element is '500[1].opbs' while the array was '500[].opbs'"""
    system = _miyabi(tmpdir, ntask=6, njob=2, pack=True)
    path = os.path.join(system._markers, "500")
    os.makedirs(path)
    with open(os.path.join(path, "chunks.json"), "w") as f:
        json.dump({"0": [0, 1, 2], "1": [3, 4, 5]}, f)
    for marker in ["0.ok", "1.ok", "2.ok", "3.fail", "4.ok", "5.ok"]:
        open(os.path.join(path, marker), "w").close()

    with patch.object(system, "monitor_job_status") as mock_monitor:
        assert(system.job_outcome("500[].opbs") == (False, "3"))
    mock_monitor.assert_not_called()


def test_pbs_entry_point_runs_the_chunk_script():
    """The script that a job starts passes what 'run_chunk' needs"""
    script = os.path.join(os.path.dirname(Chinook.run_chunk), "pbs_entry_point")
    with open(script) as f:
        text = f.read()
    assert("run_chunk)" in text)
    for variable in ["SEISFLOWS_CHUNKS", "SEISFLOWS_TASKTIME",
                     "SEISFLOWS_LOGDIR", "SEISFLOWS_MARKERS"]:
        assert(f"${{{variable}}}" in text)


def test_check_pack_is_refused_where_it_is_not_implemented(tmpdir):
    from seisflows.system.cluster import Cluster
    with open(os.path.join(str(tmpdir), "parameters.yaml"), "w") as f:
        f.write("# parameters\n")
    system = Cluster(workdir=str(tmpdir), path_system=str(tmpdir), ntask=2,
                     mpiexec="ls")
    system.pack, system.njob = True, 1
    with pytest.raises(AssertionError) as e:
        system.check()
    assert("does not support packing" in str(e.value))
