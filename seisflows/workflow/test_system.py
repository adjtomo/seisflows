#!/usr/bin/env python3
"""
This is a SeisFlows Test workflow class which is used to test out the underlying
machinery for a given set of modules, before submitting a full workflow. Used
for debugging and development, as well as ensuring that SeisFlows is in working
order.
"""
import os
from seisflows import logger
from seisflows.tools import unix
from seisflows.tools.config import Dict, get_task_id


class TestSystem:
    """
    TestSystem Workflow
    -------------
    Test individual sub-modules in a 'live' testing environment

    Parameters
    ----------

    Paths
    -----
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
        return []

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

    def run(self):
        """
        Use the system sub-module to submit some simple functions to ensure that
        we can run jobs on the system and that the job checking works
        """
        logger.info("running test workflow")

        if not "system" in self._modules:
            logger.warning("No `system` module chosen, skipping "
                           "`test_system_run`")

        system = self._modules["system"]
        system.run(funcs=[self.test_function_print_hello_world], single=False)
        system.run(funcs=[self.test_function_wait_ten_seconds], single=True)

    def test_function_print_hello_world(self):
        """Simple print function to be called by system run"""
        print(f"hello world from task id: {get_task_id()}")

    def test_function_wait_ten_seconds(self):
        """Simple wait function to be called by system run, used to test
        job status check"""
        for i in range(11):
            print(f"hello world from task id: {get_task_id()} after {i} sec.")
