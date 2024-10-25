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
from seisflows.tools.config import Dict, get_task_id


class TestFlow:
    """
    TestFlow Workflow
    -----------------
    Test individual sub-modules in a 'live' testing environment in order to
    ensure SeisFlows works appropriately given an established system and solver.

    Current functionality:
        - Submit a test function ('hello world') to the System
        - Monitor job queue for successful job completion
        - Submit an intentionally failing job and catch job failure

    .. note::

        You do not need to set System parameters `ntask`, `nproc`, `tasktime`,
        `walltime`. These will be overwritten by the setup task.

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
    ***
    """
    def __init__(self, modules=None, workdir=os.getcwd(), path_output=None,
                 **kwargs):
        """Test workflow"""
        self._modules = modules

        self.path = Dict(
            workdir=workdir,
            scratch=os.path.join(workdir, "scratch"),
            output=path_output or os.path.join(workdir, "output")
        )

    @property
    def task_list(self):
        """
        A task list which includes tests for most of the modules depending on
        whether they're included in the module list or not.
        """
        _task_list = []
        if not self.system:
            logger.warning("No `system` module chosen, skipping system tests")
        else:
            _task_list += [self.test_system_print_hello_world,
                           self.test_single_job_failure,
                           self.test_partial_array_job_failure,
                           self.test_array_job_rerun
                           ]

        return _task_list

    def check(self):    
        """
        Run check functions for all underlying modules
        """
        logger.info("running check for test workflow")
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

        # Force some internal module variables to keep testing lightweight
        logger.info("overwriting internal System parameters from given values")
        # self.system.ntask = 3
        self.system.nproc = 1
        self.system.tasktime = .25  # 15 seconds
        self.system.walltime = 2.5  # 2.5 minutes
        self.system.rerun = 0

    def run(self):
        """
        Run through the task list which should consist of various test functions
        which are meant to ensure SeisFlows works in a live working environment
        without committing a large number of resources
        """
        logger.info(msg.mjr("RUNNING TEST WORKFLOW"))
        
        for func in self.task_list:
            func()

        logger.info(msg.mjr("FINISHED TEST WORKFLOW SUCCESSFULLY"))

    def test_system_print_hello_world(self):
        """
        Use the system to run a simple print function which names the currently
        running task id to check that task id printing works. Check that the 
        output log messages show the correct task id and log statement
        """
        logger.info(msg.mnr("test array job submission"))

        def _test_function(**kwargs):
            print(f"hello world from task id: {get_task_id()}")

        # Clear the log files dir. as this is how we will check results
        unix.rm(self.system.path.log_files)
        unix.mkdir(self.system.path.log_files)

        # Run an array job
        self.system.run(funcs=[_test_function], single=False)

        time.sleep(5)  # give the system a second to catch up

        # Check the results of the run 
        log_files = glob(os.path.join(self.system.path.log_files, "*"))
        assert(len(log_files) == self.system.ntask), \
            f"number of log files does not match expected number of tasks"

        for i, fid in enumerate(sorted(log_files)):
            with open(fid, "r") as f:
                assert(f"hello world from task id: {i}" in f.read()), \
                    f"log file '{fid}' does not show correct log message"

        logger.info("finished successfully")

    def test_single_job_failure(self):
        """
        Test a single job run with built in failure to ensure that single job 
        run call modifications work and can be caught.
        """
        logger.info(msg.mnr("test for single job failure catch"))

        def _test_function(**kwargs):
            for i in range(1, 11):
                time.sleep(1)
                print(f"Task ID {get_task_id()} waited {i} sec.")

            sys.exit(-1)  # intentional job failure

        # Clear the log files dir. as this is how we will check results
        unix.rm(self.system.path.log_files)
        unix.mkdir(self.system.path.log_files)

        # Catch the system exit exception that is thrown when job fails
        try:
            self.system.run(funcs=[_test_function], single=True)
            logger.critical("submitted job was expected to fail but did not")
        except SystemExit:
            logger.info("job failure catch was successful, you can safely "
                        "ignore any error logs preceding this message")
            pass 
    
        # Check the log file for job failure
        log_files = glob(os.path.join(self.system.path.log_files, "*"))
        assert(len(log_files) == 1), f"only one log file expected"

        logger.info("finished successfully")
            
    def test_partial_array_job_failure(self):
        """
        Test that partial array job failures will cause the entire main job
        to crash. Catch the resulting error to make sure Testflow can continue
        """
        logger.info(msg.mnr("test partial array job failure"))

        def _test_function(**kwargs):
            time.sleep(10)  # need to wait for queue system to catch up
            print(f"hello world from task id: {get_task_id()}")
            if get_task_id() in [1, 2]:  # assuming `ntask`==3
                sys.exit(-1)  # intentional job failure

        # Clear the log files dir. as this is how we will check results
        unix.rm(self.system.path.log_files)
        unix.mkdir(self.system.path.log_files)

        # Catch the system exit exception that is thrown when job fails
        try:
            self.system.run(funcs=[_test_function])
        except SystemExit:
            logger.info("job failure catch was successful, you can safely "
                        "ignore any error logs preceding this message")
            pass 
    
        # Check the log file for job failure
        log_files = glob(os.path.join(self.system.path.log_files, "*"))
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
        logger.info(msg.mnr("test system job failure rerun"))

        self.system.rerun = 1
        self.system.ntask = 3
        failed_task_ids = [1, 2]  # 0 should pass, 1 and 2 should fail

        def _test_function(**kwargs):
            time.sleep(10)  # need to wait for queue system to catch up
            print(f"hello world from task id: {get_task_id()}")
            if get_task_id() in failed_task_ids: 
                sys.exit(-1)  # intentional job failure

        # Clear the log files dir. as this is how we will check results
        unix.rm(self.system.path.log_files)
        unix.mkdir(self.system.path.log_files)

        # Catch the system exit exception that is thrown when job fails
        try:
            self.system.run(funcs=[_test_function])
        except SystemExit:
            logger.info("job failure catch was successful, you can safely "
                        "ignore any error logs preceding this message")
            pass 
    
        # Check the log file for job failure
        log_files = glob(os.path.join(self.system.path.log_files, "*"))
        nlog_files = self.system.ntask + len(failed_task_ids)
        assert(len(log_files) == nlog_files), f"{nlog_files} log files expected"

        logger.info("finished successfully")
