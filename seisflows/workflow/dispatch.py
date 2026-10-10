#!/usr/bin/env python3
"""
Submission of workflow tasks to the system module, blocking or not.

Everything that decides HOW a task-list function's work reaches the
scheduler lives here, so that the workflow classes can stay about the
science: a task calls `self.dispatcher.run(funcs, ...)` where it would
otherwise have called `self.system.run(funcs, ...)`.

With `system.blocking` True (default) that IS `system.run()`: the master job
waits for the work to finish. With it False the master job instead dispatches
the work, schedules a dependent job that resumes the workflow once it is
done, and exits, so it never has to outlive a single scheduler job. See
`Dispatcher.run` for the details.
"""
import os
import sys
import logging
from contextlib import contextmanager

from seisflows import logger
from seisflows.tools import msg
from seisflows.tools.config import atomic_write_lines


class WorkflowSuspended(Exception):
    """
    Raised by `Dispatcher.run()` after it dispatches work non-blockingly, to
    stop the calling task-list function (nothing after the dispatch runs).
    `Forward.run()` catches it and exits the current process cleanly: it is
    NOT a task failure, the task is simply still in progress elsewhere. A
    dependent job later re-invokes ``seisflows submit`` to resume.
    """


def pending_state(job_id, attempt=0):
    """
    Build the state-file value for a task that has dispatched work
    non-blockingly and is waiting on it

    :type job_id: str
    :param job_id: ID of the dispatched job, from `system.run(blocking=False)`
    :type attempt: int
    :param attempt: how many times the failed tasks of this job have already
        been re-dispatched (0 for the original dispatch), used to honor
        `system.rerun` across process restarts
    :rtype: str
    :return: state value of the form 'PENDING:<job_id>:<attempt>'
    """
    return f"PENDING:{job_id}:{attempt}"


def parse_pending_state(val):
    """
    Inverse of `pending_state`. Tolerates a missing attempt count (e.g.
    'PENDING:12345') as attempt 0.

    :type val: str
    :param val: state-file value beginning with 'PENDING'
    :rtype: (str, int)
    :return: (job_id, attempt)
    :raises ValueError: if `val` is not a well-formed PENDING value
    """
    tag, job_id, *rest = val.split(":")
    if tag != "PENDING" or not job_id or len(rest) > 1:
        raise ValueError(f"malformed PENDING state: '{val}'")
    return job_id, int(rest[0]) if rest else 0


class ResumeFilter(logging.Filter):
    """
    Logging filter that lets only warnings and worse through. Installed while
    a resumed workflow replays its startup, so the log of the whole workflow
    reads as one continuous story rather than repeating it at every resume.
    """
    def filter(self, record):
        """
        :type record: logging.LogRecord
        :param record: the record being logged
        :rtype: bool
        :return: True if the record should be logged
        """
        return record.levelno >= logging.WARNING


class Dispatcher:
    """
    Runs a workflow's tasks on its system module, deciding between blocking
    and non-blocking submission.

    .. note::
        Reads and writes each task's state through the workflow it belongs to
        (`workflow._states`, `workflow.checkpoint()`) rather than owning any
        state itself: the state file records the progress of every task,
        blocking or not, so it stays with the workflow.

    :type workflow: seisflows.workflow.forward.Forward
    :param workflow: the workflow whose tasks are being dispatched. Needs
        `system`, `path`, `_states` and `checkpoint()`
    """
    def __init__(self, workflow):
        self.workflow = workflow
        self._task = None  # name of the task-list function now running
        self._key = None  # key set by `keyed()`, for dispatches inside it
        self._counts = {}  # dispatches made so far by that task, per key

    @contextmanager
    def task(self, name):
        """
        Context in which the task-list function `name` runs, so that
        dispatches made inside it are attributed to that task. Used by
        `Forward.run()`.

        :type name: str
        :param name: name of the task-list function
        """
        self._task, self._counts = name, {}
        try:
            yield
        finally:
            self._task = None

    @contextmanager
    def keyed(self, name):
        """
        Context in which dispatches are recorded under the key `name`
        instead of the running task's name. For tasks that reach their
        dispatch through other methods (e.g. a loop calling
        `super().evaluate_initial_misfit()` once per force), where passing a
        `key` argument down would mean changing every signature in between.

        :type name: str
        :param name: state-file key for dispatches made inside the context
        """
        self._key = name
        try:
            yield
        finally:
            self._key = None

    @property
    def nonblocking(self):
        """
        True if work dispatched right now will be non-blocking: a
        task-list function is running (so there is something to suspend and
        resume) and the system module is set to `blocking=False`.

        :rtype: bool
        :return: whether dispatches are currently non-blocking
        """
        return (self._task is not None and
                not getattr(self.workflow.system, "blocking", True))

    @property
    def resuming(self):
        """
        True if the running task-list function already dispatched work
        non-blockingly in an earlier process and is being re-entered to
        collect it. Anything a task does BEFORE dispatching that is not safe
        to repeat (e.g. deleting directories the dispatched work writes
        into) must be skipped when this is True.

        :rtype: bool
        :return: whether this task is resuming a non-blocking dispatch
        """
        return self._pending(self._task)

    def begin_resume(self):
        """
        Called once the workflow has read its state file. If a dispatch is
        waiting to be collected, this process is a resumption of an earlier
        one, and everything it logs until it gets back to that dispatch
        (startup, banners of finished tasks, whatever the task does before
        dispatching) repeats what the log already says. Hide it, other than
        warnings and errors. `end_resume` shows the log again.
        """
        self.end_resume()
        if any(self._pending(key) for key in self.workflow._states):
            logger.addFilter(ResumeFilter())

    def end_resume(self, key=None, job_id=None):
        """
        Show the log again after `begin_resume`. Does nothing if that did
        not hide it. Always safe to call.

        :type key: str
        :param key: state-file key of the dispatch being collected, if the
            log is being shown again because the workflow got back to one
        :type job_id: str
        :param job_id: ID of the job that was dispatched under `key`
        """
        hidden = [f for f in logger.filters if isinstance(f, ResumeFilter)]
        for resume_filter in hidden:
            logger.removeFilter(resume_filter)

        if hidden and key is not None:
            iteration = getattr(self.workflow, "iteration", None)
            where = f" ITERATION {iteration:0>2}" if iteration else ""
            logger.info(f"---- RESUMING{where}: '{self._task}' after job "
                        f"{job_id} ----")
            logger.debug(f"collecting dispatch '{key}'")

    def _pending(self, key):
        """
        Whether the state recorded under `key` is a PENDING dispatch

        :type key: str
        :param key: state-file key of a dispatch, see `run`
        :rtype: bool
        :return: True if a dispatch under `key` is waiting to be collected
        """
        state = self.workflow._states.get(key)
        return isinstance(state, str) and state.startswith("PENDING")

    @property
    def marker(self):
        """
        Path to a file that a failed compute-side preflight check writes its
        error into (see `flag_preflight_failure`). One file shared by all
        tasks, since they all check the same thing: those that fail write
        identical content. Lives in scratch, visible to master and compute
        nodes alike.

        :rtype: str
        :return: path to the marker file
        """
        return os.path.join(self.workflow.path.scratch, "preflight_failed.txt")

    def run(self, funcs, key=None, **kwargs):
        """
        Run `funcs` on the system module, non-blockingly if
        `system.blocking` is False. Task-list functions call this instead of
        `system.run()`.

        Blocking (default): identical to `system.run(funcs, **kwargs)`.

        Non-blocking, the work happens in two process lifetimes, and the
        calling task-list function is simply entered twice:

        1. First entry: dispatch `funcs`, record the job ID as the 'PENDING'
           state of this task (or of `key`), schedule a dependent job that
           resumes the workflow once the dispatched job finishes, and
           suspend (raises `WorkflowSuspended`, so nothing after this call
           runs).
        2. Re-entry in the resumed process: the dispatched job has
           finished. Check its outcome. On success, return so the calling
           function can carry on with whatever follows the call (reading
           results, etc.). On failure, re-dispatch only the failed tasks
           (up to `system.rerun` times, like blocking `system.run()` does)
           and suspend again, or exit if retries are exhausted. A failure
           flagged by a preflight check is never retried, see
           `flag_preflight_failure`.

        .. warning::
            Because the caller is re-entered, everything it does before
            this call runs again on resume and must be safe to repeat (see
            `resuming` for the exceptions).

        .. note::
            Each dispatch's PENDING state is kept under a key in the state
            file, by default the task's own name. A task that dispatches
            several times (e.g. once per line search step) must give each
            dispatch its own `key`; two dispatches under one key in a
            single invocation are refused, since they could not be told
            apart on resume. Keys must not contain ':' (the state file
            splits each line at its first one). The task is responsible for
            resetting a key it wants to reuse, e.g. after a restart.

        :type funcs: list of methods
        :param funcs: functions to run on the system, see `system.run()`
        :type key: str
        :param key: state-file key to record this dispatch under. Defaults to
            the key of the enclosing `keyed()` context if there is one, else
            the name of the running task-list function. Not passed on to
            `system.run()`
        :param kwargs: passed through to `system.run()`
        :raises WorkflowSuspended: after dispatching, see above
        """
        system = self.workflow.system
        if not self.nonblocking:
            return system.run(funcs, **kwargs)

        key = key or self._key or self._task
        self._counts[key] = self._counts.get(key, 0) + 1
        assert(self._counts[key] == 1), (
            f"task '{self._task}' dispatched more than once under the key "
            f"'{key}' while `system.blocking` is False, which is not "
            f"supported: its PENDING state cannot distinguish between "
            f"multiple dispatches. Give each dispatch its own `key`"
        )

        if not self._pending(key):
            self.end_resume()
            # A marker left by an earlier, since-fixed failure must not be
            # mistaken for a failure of this dispatch
            self.clear_preflight_failure()
            job_id = system.run(funcs, blocking=False, **kwargs)
            self._suspend(key, job_id, attempt=0)

        job_id, attempt = parse_pending_state(self.workflow._states[key])
        self.end_resume(key, job_id)
        success, failed = system.job_outcome(job_id)
        if success:
            logger.info(f"job {job_id} finished successfully")
            return

        # A failed preflight check is deterministic, retrying would only
        # fail the same way, so it ends the workflow instead of using up
        # `system.rerun` attempts
        self._exit_on_preflight_failure(key, job_id)

        if attempt < system.rerun:
            logger.warning(f"job {job_id} had failed tasks: {failed}. "
                           f"Attempt {attempt + 1}/{system.rerun} to rerun")
            job_id = system.run(funcs, blocking=False, array=failed, **kwargs)
            self._suspend(key, job_id, attempt=attempt + 1)

        # Out of retries. Mark failed (not PENDING) so that fixing the
        # problem and re-running `seisflows submit` dispatches fresh, as
        # happens after a failure in blocking mode
        self.workflow._states[key] = -1
        self.workflow.checkpoint()
        logger.critical(
            msg.cli(f"Stopping workflow. Please check logs for details",
                    items=[f"TASK:   {key}", f"JOB_ID: {job_id}",
                           f"FAILED TASK IDS: {failed}"],
                    header="system run error", border="=")
        )
        sys.exit(-1)

    def flag_preflight_failure(self, error):
        """
        Record that a check run on a compute node failed. Called by the
        check itself, just before it re-raises to fail its task. The master
        job finds this after the job fails and, instead of retrying a
        failure that would only repeat, exits with the check's message.

        :type error: Exception
        :param error: the exception the check raised
        """
        atomic_write_lines(self.marker, [f"{type(error).__name__}: {error}\n"])

    def clear_preflight_failure(self):
        """Remove a preflight failure marker, if there is one"""
        if os.path.exists(self.marker):
            os.remove(self.marker)

    def _exit_on_preflight_failure(self, key, job_id):
        """
        If a compute-side preflight check failed (see
        `flag_preflight_failure`), mark `key` failed and exit with the
        check's own error message. Does nothing if there is no failure
        marker, i.e. the job failed for some other (possibly transient)
        reason.

        :type key: str
        :param key: state-file key of the dispatch whose job failed
        :type job_id: str
        :param job_id: ID of the job that failed
        """
        if not os.path.exists(self.marker):
            return

        with open(self.marker, "r") as f:
            error = f.read().strip()
        self.clear_preflight_failure()

        # Failed, not PENDING: once the cause is fixed, `seisflows submit`
        # then dispatches fresh
        self.workflow._states[key] = -1
        self.workflow.checkpoint()
        logger.critical(
            msg.cli(f"A check run on the compute nodes failed. The job is "
                    f"not retried because it would fail the same way.",
                    items=[f"TASK:   {key}", f"JOB_ID: {job_id}",
                           f"ERROR:  {error}"],
                    header="preflight check error", border="=")
        )
        sys.exit(-1)

    def _suspend(self, key, job_id, attempt):
        """
        Record that the dispatch under `key` is waiting on `job_id`, schedule
        the dependent job that resumes the workflow when it finishes, and
        suspend.

        :type key: str
        :param key: state-file key of the dispatch that created `job_id`
        :type job_id: str
        :param job_id: ID of the job the workflow is waiting on
        :type attempt: int
        :param attempt: number of re-dispatches so far, see `pending_state`
        :raises WorkflowSuspended: always
        """
        # Checkpoint BEFORE scheduling the resume job so that it can never
        # start and read a state file that doesn't yet know about `job_id`
        self.workflow._states[key] = pending_state(job_id, attempt)
        self.workflow.checkpoint()
        resume_id = self.workflow.system.submit(
            workdir=self.workflow.path.workdir,
            parameter_file=self.workflow.system.path.par_file,
            depends_on=job_id
        )

        # The log of the master job is the high-level story, point to where
        # the output of the dispatched tasks goes
        log_files = self.workflow.system.path.log_files
        if log_files.startswith(self.workflow.path.workdir):
            log_files = os.path.relpath(log_files, self.workflow.path.workdir)
        logger.info(f"dispatched job {job_id}, job {resume_id} will resume "
                    f"the workflow when it finishes. Task logs are in "
                    f"{log_files}")
        raise WorkflowSuspended()
