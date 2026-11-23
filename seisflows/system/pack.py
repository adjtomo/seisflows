#!/usr/bin/env python3
"""
Pack mode: run the tasks of an array in fewer jobs.

A scheduler treats every element of an array as a job of its own, and each one
can wait in the queue. When the queue is long (or only a few jobs may run at
once) it is quicker to ask for a few jobs that each run their share of the
tasks one after the other. The system modules do this when `pack` is True,
splitting the tasks into `njob` blocks of neighbouring ids, one block per
array element, run by the `run_chunk` script.

Each task is run in a process of its own by the usual `run` script, so one
task cannot leak its working directory, environment or exit into the next, and
it has its own log and time limit. A task leaves a marker file when it is
done, `<markers>/<job id>/<task id>.ok` or `.fail`, which is how a job (and
any task that never started because the job ended early) is known to have
failed, rather than the state that the scheduler reports for the whole block.
"""
import json
import math
import os
import re
import subprocess
import sys

from seisflows.tools.config import get_task_id

# Extra time requested over the sum of the task times, for the parts of a job
# that are not tasks (starting the job, loading the environment)
PACK_BUFFER = 0.1

# The script that runs a single task
RUN_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "runscripts", "run")


def parse_array(array):
    """
    Task ids from a scheduler array string, e.g., '0,2,5-7' or '0-79%15'. A
    throttle ('%15') is ignored, as it only limits how many run at once.

    :type array: str
    :param array: array string, comma separated ids and ranges 'start-end'
    :rtype: list of int
    :return: task ids in the order given
    """
    task_ids = []
    for item in str(array).split("%")[0].split(","):
        item = item.strip()
        if not item:
            continue
        start, _, end = item.partition("-")
        task_ids.extend(range(int(start), int(end or start) + 1))
    return task_ids


def split_tasks(task_ids, njob):
    """
    Split task ids into at most `njob` blocks of neighbouring ids, as even as
    they can be, e.g., 10 tasks in 4 blocks are 3, 3, 2 and 2. The block
    sizes of tasks that cost the same differ by at most one task.

    :type task_ids: list of int
    :param task_ids: task ids to split, in order
    :type njob: int
    :param njob: number of blocks. If it is more than the number of tasks,
        every task is its own block
    :rtype: list of list of int
    :return: the blocks, the largest first
    """
    njob = max(1, min(int(njob), len(task_ids)))
    size, extra = divmod(len(task_ids), njob)
    blocks, start = [], 0
    for i in range(njob):
        end = start + size + (1 if i < extra else 0)
        blocks.append(list(task_ids[start:end]))
        start = end
    return blocks


def pack_tasktime(tasktime, ntask):
    """
    Time limit of a job that runs `ntask` tasks one after the other

    :type tasktime: float
    :param tasktime: time limit of one task in minutes
    :type ntask: int
    :param ntask: number of tasks the job runs, i.e., the size of the largest
        block
    :rtype: int
    :return: time limit in whole minutes
    """
    # Rounded first, 10 x 18 x 1.1 is 198.00000000000003
    return int(math.ceil(round(tasktime * ntask * (1 + PACK_BUFFER), 6)))


def job_key(job_id):
    """
    The part of a job id that is the same for every element of an array and
    for the job id that the scheduler returned when the array was submitted,
    e.g., '123' for '123[].opbs', '123[4].opbs' and '123_4'

    :type job_id: str
    :param job_id: a job id, from the scheduler or from the job's environment
    :rtype: str
    :return: leading digits of the job id, the whole id if there are none
    """
    match = re.match(r"\d+", str(job_id).strip())
    return match.group(0) if match else str(job_id).strip()


def failed_tasks(markers, job_id):
    """
    Which tasks of a packed job did not finish. Tasks that left no marker
    count as failed, the job may have ended before it got to them.

    :type markers: str
    :param markers: directory of the markers of all packed jobs
    :type job_id: str
    :param job_id: id of the packed job, as returned when it was submitted
    :rtype: list of int or NoneType
    :return: ids of the tasks that did not finish, None if the job was not
        packed (it has no list of tasks)
    """
    path = os.path.join(markers, job_key(job_id))
    if not os.path.exists(os.path.join(path, "chunks.json")):
        return None
    with open(os.path.join(path, "chunks.json")) as f:
        chunks = json.load(f)
    return sorted(task for chunk in chunks.values() for task in chunk
                  if not os.path.exists(os.path.join(path, f"{task}.ok")))


def run_chunk(chunk, funcs, kwargs, tasktime, logdir, markers, job_id,
              environment=None, run_script=RUN_SCRIPT):
    """
    Run the tasks of one block one after the other, each by the `run` script
    in a process of its own with the task id of the task. A task that fails or
    runs for longer than `tasktime` does not stop the ones after it.

    :type chunk: list of int
    :param chunk: ids of the tasks to run
    :type funcs: str
    :param funcs: path to the pickle file of functions, as given to `run`
    :type kwargs: str
    :param kwargs: path to the pickle file of keyword arguments, as given to
        `run`
    :type tasktime: float
    :param tasktime: time limit of one task in minutes
    :type logdir: str
    :param logdir: directory for the log of each task, `<job>_task<id>.log`
    :type markers: str
    :param markers: directory of the markers of all packed jobs
    :type job_id: str
    :param job_id: id of the packed job that this is an element of
    :type environment: str
    :param environment: environment variables for the tasks, as given to `run`
    :type run_script: str
    :param run_script: the script that runs one task
    :rtype: int
    :return: 0 if every task finished, 1 if any did not
    """
    key = job_key(job_id)
    path = os.path.join(markers, key)
    os.makedirs(path, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    command = [sys.executable, run_script, "--funcs", funcs, "--kwargs",
               kwargs]
    if environment:
        command += ["--environment", environment]

    exit_code = 0
    for task in chunk:
        # `SEISFLOWS_TASKID` is read before the array index of the scheduler
        env = dict(os.environ, SEISFLOWS_TASKID=str(task))
        with open(os.path.join(logdir, f"{key}_task{task}.log"), "w") as log:
            try:
                result = subprocess.run(command, env=env, stdout=log,
                                        stderr=subprocess.STDOUT,
                                        timeout=tasktime * 60)
                done = result.returncode == 0
                if not done:
                    log.write(f"\ntask {task} exited with code "
                              f"{result.returncode}\n")
            except subprocess.TimeoutExpired:
                done = False
                log.write(f"\ntask {task} exceeded its time limit of "
                          f"{tasktime} minutes\n")
        open(os.path.join(path, f"{task}.{'ok' if done else 'fail'}"),
             "w").close()
        exit_code = exit_code or int(not done)

    return exit_code


def main():
    """
    Entry point of the `run_chunk` script: run the block of tasks that the
    array element (its index, see `get_task_id`) is responsible for
    """
    import argparse

    parser = argparse.ArgumentParser("Run a block of tasks of a packed job")
    parser.add_argument("--funcs", required=True)
    parser.add_argument("--kwargs", required=True)
    parser.add_argument("--chunks", required=True,
                        help="path to the JSON file of blocks of task ids")
    parser.add_argument("--tasktime", type=float, required=True,
                        help="time limit of one task in minutes")
    parser.add_argument("--logdir", required=True)
    parser.add_argument("--markers", required=True)
    parser.add_argument("--environment", nargs="?", default=None)
    args = parser.parse_args()

    with open(args.chunks) as f:
        chunk = json.load(f)[str(get_task_id())]
    job_id = (os.getenv("SLURM_ARRAY_JOB_ID") or os.getenv("PBS_JOBID") or
              "local")
    sys.exit(run_chunk(chunk, args.funcs, args.kwargs, args.tasktime,
                       args.logdir, args.markers, job_id, args.environment))
