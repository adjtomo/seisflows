#!/usr/bin/env python3
"""
This is a SeisFlows Test workflow class which is used to test out the underlying
machinery for a given set of modules, before submitting a full workflow. Used
for debugging and development, as well as ensuring that SeisFlows is in working
order before committing a large number of compute resources.
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
            _task_list.append(self.test_system_print_hello_world)
            _task_list.append(self.test_system_wait_ten_seconds_then_fail)

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
        self.system.ntask = 3
        self.system.nproc = 1
        self.system.tasktime = .25  # 15 seconds
        self.system.walltime = 2.5  # 2.5 minutes

    def run(self):
        """
        Run through the task list which should consist of various test functions
        which are meant to ensure SeisFlows works in a live working environment
        without committing a large number of resources
        """
        logger.info(msg.mjr("RUNNING TEST WORKFLOW"))
        
        for func in self.task_list:
            func()

        logger.info(msg.mjr("FINISHED TEST WORKFLOW"))

    def test_system_print_hello_world(self):
        """
        Use the system to run a simple print function which names the currently
        running task id to check that task id printing works. Check that the 
        output log messages show the correct task id and log statement
        """
        logger.info("running system test for job array submission")

        def _test_function():
            print(f"hello world from task id: {get_task_id()}")

        # Clear the log files dir. as this is how we will check results
        unix.rm(self.system.path.log_files)
        unix.mkdir(self.system.path.log_files)

        # Run an array job
        self.system.run(funcs=[_test_function], single=False)

        time.sleep(5)  # give the system a second to catch up

        # Check the results of the run 
        log_files = glob(os.path.join(self.system.path.log_files, "*_*"))
        assert(len(log_files) == self.system.ntask), \
            f"number of log files does not match expected number of tasks"

        for i, fid in enumerate(sorted(log_files)):
            with open(fid, "r") as f:
                assert(f"hello world from task id: {i}" in f.read()), \
                    f"log file '{fid}' does not show correct log message"

        logger.info("job array submission system test finished successfully")
                

    def test_system_wait_ten_seconds_then_fail(self):
        """
        Simple wait function to be called by system run, used to test
        job status check function which monitors job queue. Throw in a job
        failure at the end of the function to check that system exit works
        """
        logger.info("running system test for job queue monitoring and "
                    "job failure catching")

        def _test_function():
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
            logger.critical("job was expected to fail but did not")
        except SystemExit:
            logger.info("job failure catch was successful")
            pass 
    
        # Check the log file for job failure
        log_files = glob(os.path.join(self.system.path.log_files, "*_*"))
        assert(len(log_files) == 1), f"only one log file expected"

        logger.info("job queue and fail system test finished successfully")
            
