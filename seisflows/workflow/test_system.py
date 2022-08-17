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


class Test:
    """
    Test Workflow
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

    def setup(self):
        """
        Creates required directory structure
        """
        for path in [self.path.workdir, self.path.scratch, self.path.output]:
            unix.mkdir(path)

        for module in self._modules:
            module.setup()

    def test_system_run(self):
        """
        Use the system sub-module to submit some simple functions to ensure that
        we can run jobs on the system and that the job checking works
        """
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
