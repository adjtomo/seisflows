#!/usr/bin/env python3
"""
This is a SeisFlows Test workflow class which is used to test out the underlying
`System` machinery for a given set of modules.

TestFlow is intended for development and debugging purposes. It is used to
ensure that these tasks are not done on large-scale, compute-intensive jobs.
"""
import os
import time
import sys
from glob import glob
from seisflows import logger
from seisflows.tools import unix, msg
from seisflows.tools.config import get_task_id, atomic_write_lines
from seisflows.workflow.dispatch import parse_pending_state
from seisflows.workflow.forward import Forward


class TestFlow(Forward):
    """
    TestFlow Workflow
    -----------------
    Test individual sub-modules in a 'live' testing environment in order to
    ensure SeisFlows works appropriately given an established system and solver.

    Current functionality:
        - Submit a test function ('hello world') to the System
        - Monitor job queue for successful job completion
        - Submit an intentionally failing job and catch job failure
        - Re-run the failed tasks of an array job
        - If `system.blocking` is False: check that the main job is resumed
          by the scheduler after each of these (so, the scheduler's job
          dependencies), that the queue and run times of jobs are available,
          and that a failed check on a compute node is not retried

    .. note::

        You do not need to set System parameters `ntask`, `nproc`, `tasktime`,
        `walltime`. These will be overwritten by the setup task.

    .. note::

        Nothing needs to be compiled: the solver, preprocessing and
        optimization modules are not used, so they can be 'null'. With
        `system.blocking` False each test is a trip through the queue and
        the main job is restarted after each, so it takes longer to finish.

    Parameters
    ----------

    Paths
    -----
    :type workdir: str
    :param workdir: working directory in which to perform a SeisFlows workflow.
        SeisFlows internal directory structure will be created here. Default cwd
    :type path_output: str
    :param path_output: path to directory used for permanent storage on disk.
        Results and expored scratch files are saved here.
    :type path_state_file: str
    :param path_state_file: path to a text file used to track the current
        status of a workflow (i.e., what functions have already been completed),
        used for resuming the workflow after a non-blocking job
    ***
    """
    # Not a test for pytest to collect, despite the name
    __test__ = False

    def __init__(self, modules=None, workdir=os.getcwd(), path_output=None,
                 path_state_file=None, **kwargs):
        """Test workflow"""
        super().__init__(modules=modules, workdir=workdir,
                         path_output=path_output,
                         path_state_file=path_state_file)

    @property
    def task_list(self):
        """
        A task list which includes tests for most of the modules depending on
        whether they're included in the module list or not.

        .. note::
            The parent class reads this when it is made, before `setup()`
            assigns the modules, so it looks in the list of modules (which is
            also None when a parameter file is made).
        """
        system = self._modules["system"] if self._modules else None

        _task_list = []
        if system:
            _task_list += [self.test_system_print_hello_world,
                           self.test_single_job_failure,
                           self.test_partial_array_job_failure,
                           self.test_array_job_rerun
                           ]
            # Only non-blocking jobs look for a check that failed on a node
            if not getattr(system, "blocking", True):
                _task_list.append(self.test_preflight_failure_is_not_retried)

        return _task_list

    def check(self):
        """
        Run check functions for all underlying modules
        """
        logger.info("running check for test workflow")
        if not self._modules["system"]:
            logger.warning("No `system` module chosen, skipping system tests")
        for name, module in self._modules.items():
            if module:
                module.check()

    def setup(self):
        """
        Creates required directory structure
        """
        logger.info("running setup for test workflow")
        for path in [self.path.workdir, self.path.scratch, self.path.output]:
            unix.mkdir(path)

        for name, module in self._modules.items():
            if module:
                module.setup()

        # Assign modules to internal attributes, some may be Null
        self.system = self._modules["system"]
        self.solver = self._modules["solver"]
        self.preprocess = self._modules["preprocess"]
        self.optimize = self._modules["optimize"]

        # Keeps track of which tests have been completed, which a non-blocking
        # main job needs to carry on from where it was
        if not os.path.exists(self.path.state_file):
            with open(self.path.state_file, "w") as f:
                f.write("# SeisFlows State File\n")

        # Force some internal module variables to keep testing lightweight
        logger.info("overwriting internal System parameters from given values")
        self.system.ntask = 3
        self.system.nproc = 1
        self.system.tasktime = 1  # 15 seconds
        self.system.walltime = 5  # 2.5 minutes
        self.system.rerun = 0

    def run(self):
        """
        Run through the task list which should consist of various test functions
        which are meant to ensure SeisFlows works in a live working environment
        without committing a large number of resources

        .. note::
            When every test has passed the state file is reset, because unlike
            other workflows the point of this one is to be run again, e.g.,
            after changing the parameters of the system.

        :rtype: str
        :return: 'COMPLETE' if every test passed, 'SUSPENDED' if one is waiting
            for a job, see `workflow.forward.Forward.run`
        """
        logger.info(msg.mjr("RUNNING TEST WORKFLOW"))

        status = super().run()

        if status == "COMPLETE":
            logger.info(msg.mjr("FINISHED TEST WORKFLOW SUCCESSFULLY"))
            self._states = {key: 0 for key in self._states}
            self.checkpoint()

        return status

    def test_system_print_hello_world(self):
        """
        Use the system to run a simple print function which names the currently
        running task id to check that task id printing works. Check that the
        output log messages show the correct task id and log statement. If the
        main job is non-blocking, also check that the times that the job
        queued and ran for are available.
        """
        logger.info("test array job submission")

        def _test_function(**kwargs):
            print(f"hello world from task id: {get_task_id()}")

        # Clear the log files dir. as this is how we will check results. A
        # non-blocking job comes back to this function after the job is run,
        # so only do this the first time
        if not self._dispatcher.resuming:
            unix.rm(self.system.path.log_files)
            unix.mkdir(self.system.path.log_files)

        # Run an array job
        self._dispatcher.run(funcs=[_test_function], single=False)

        time.sleep(5)  # give the system a second to catch up

        # Check the results of the run, other than copies of the main log
        log_files = [fid for fid in glob(os.path.join(
                         self.system.path.log_files, "*")) if not
                     os.path.basename(fid).startswith(("sflog", "parameters"))]
        assert(len(log_files) == self.system.ntask), \
            f"number of log files does not match expected number of tasks"

        for i, fid in enumerate(sorted(log_files)):
            with open(fid, "r") as f:
                assert(f"hello world from task id: {i}" in f.read()), \
                    f"log file '{fid}' does not show correct log message"

        # The job's ID is how a non-blocking job remembers what it waits for
        state = self._states.get("test_system_print_hello_world")
        if isinstance(state, str) and state.startswith("PENDING"):
            job_id, _ = parse_pending_state(state)
            times = self.system.job_times(job_id)
            if times is None:
                logger.warning(f"system cannot say how long job {job_id} "
                               f"queued and ran for, `system.job_times()` "
                               f"does not work for this scheduler. This only "
                               f"affects the log")
            else:
                logger.info(f"job {job_id} queued for "
                            f"{msg.duration(times[0])} and ran for "
                            f"{msg.duration(times[1])}")

        logger.info("finished successfully")

    def test_single_job_failure(self):
        """
        Test a single job run with built in failure to ensure that single job
        run call modifications work and can be caught.
        """
        logger.info("test for single job failure catch")

        def _test_function(**kwargs):
            for i in range(1, 11):
                time.sleep(1)
                print(f"Task ID {get_task_id()} waited {i} sec.")

            sys.exit(-1)  # intentional job failure

        # Clear the log files dir. as this is how we will check results
        if not self._dispatcher.resuming:
            unix.rm(self.system.path.log_files)
            unix.mkdir(self.system.path.log_files)

        # Catch the system exit exception that is thrown when job fails
        failed = False
        try:
            self._dispatcher.run(funcs=[_test_function], single=True)
        except SystemExit:
            failed = True
            logger.info("job failure catch was successful, you can safely "
                        "ignore any error logs preceding this message")
        assert(failed), "submitted job was expected to fail but did not"

        # Check the log file for job failure
        log_files = [fid for fid in glob(os.path.join(
                         self.system.path.log_files, "*")) if not
                     os.path.basename(fid).startswith(("sflog", "parameters"))]
        assert(len(log_files) == 1), f"only one log file expected"

        logger.info("finished successfully")

    def test_partial_array_job_failure(self):
        """
        Test that partial array job failures will cause the entire main job
        to crash. Catch the resulting error to make sure Testflow can continue
        """
        logger.info("test partial array job failure")

        def _test_function(**kwargs):
            time.sleep(10)  # need to wait for queue system to catch up
            print(f"hello world from task id: {get_task_id()}")
            if get_task_id() in [1, 2]:  # assuming `ntask`==3
                sys.exit(-1)  # intentional job failure

        # Clear the log files dir. as this is how we will check results
        if not self._dispatcher.resuming:
            unix.rm(self.system.path.log_files)
            unix.mkdir(self.system.path.log_files)

        # Catch the system exit exception that is thrown when job fails
        failed = False
        try:
            self._dispatcher.run(funcs=[_test_function])
        except SystemExit:
            failed = True
            logger.info("job failure catch was successful, you can safely "
                        "ignore any error logs preceding this message")
        assert(failed), "submitted job was expected to fail but did not"

        # Check the log file for job failure
        log_files = [fid for fid in glob(os.path.join(
                         self.system.path.log_files, "*")) if not
                     os.path.basename(fid).startswith(("sflog", "parameters"))]
        assert(len(log_files) == 3), f"unexpected number of log files found"

        logger.info("finished successfully")

    def test_array_job_rerun(self):
        """
        Test that partial array job failures will cause the entire main job
        to crash. The idea behind this function is that we want a job to fail
        once, and then pass on the second attempt. Because of the architecture
        of the failure recovery mechanism, we cannot incremenet a counter to
        check if the job has been rereun, so we check that the correct amount
        of log files has been produced
        """
        logger.info("test system job failure rerun")

        # Set again every time as a non-blocking job runs setup() each time
        self.system.rerun = 1
        self.system.ntask = 3
        failed_task_ids = [1, 2]  # 0 should pass, 1 and 2 should fail

        def _test_function(**kwargs):
            time.sleep(10)  # need to wait for queue system to catch up
            print(f"hello world from task id: {get_task_id()}")
            if get_task_id() in failed_task_ids:
                sys.exit(-1)  # intentional job failure

        # Clear the log files dir. as this is how we will check results
        if not self._dispatcher.resuming:
            unix.rm(self.system.path.log_files)
            unix.mkdir(self.system.path.log_files)

        # Catch the system exit exception that is thrown when job fails
        failed = False
        try:
            self._dispatcher.run(funcs=[_test_function])
        except SystemExit:
            failed = True
            logger.info("job failure catch was successful, you can safely "
                        "ignore any error logs preceding this message")
        assert(failed), "submitted job was expected to fail but did not"

        # Check the log file for job failure
        log_files = [fid for fid in glob(os.path.join(
                         self.system.path.log_files, "*")) if not
                     os.path.basename(fid).startswith(("sflog", "parameters"))]
        nlog_files = self.system.ntask + len(failed_task_ids)
        assert(len(log_files) == nlog_files), f"{nlog_files} log files expected"

        logger.info("finished successfully")

    def test_preflight_failure_is_not_retried(self):
        """
        Test that a check that fails on a compute node ends the workflow
        without the failed job being run again, however many `rerun` there are.
        The function tells the main job by writing a file to the scratch
        directory, which also checks that the main job and the compute nodes
        share a file system. Only for non-blocking jobs, which are the only ones
        that look for such a file.
        """
        logger.info("test that a failed check on a compute node is not retried")

        self.system.rerun = 2  # would be 3 jobs if this was not special
        marker = self._dispatcher.marker

        def _test_function(**kwargs):
            atomic_write_lines(marker, ["RuntimeError: intentional failure\n"])
            sys.exit(-1)  # intentional job failure

        # Clear the log files dir. as this is how we will check results
        if not self._dispatcher.resuming:
            unix.rm(self.system.path.log_files)
            unix.mkdir(self.system.path.log_files)

        failed = False
        try:
            self._dispatcher.run(funcs=[_test_function], single=True)
        except SystemExit:
            failed = True
            logger.info("job failure catch was successful, you can safely "
                        "ignore any error logs preceding this message")
        assert(failed), "submitted job was expected to fail but did not"

        # One job, not the original and two re-runs, and the file was found
        log_files = [fid for fid in glob(os.path.join(
                         self.system.path.log_files, "*")) if not
                     os.path.basename(fid).startswith(("sflog", "parameters"))]
        assert(len(log_files) == 1), \
            f"only one log file expected, the failed job was run again"
        assert(not os.path.exists(marker)), \
            f"main job did not find the file written by the compute node"

        logger.info("finished successfully")
