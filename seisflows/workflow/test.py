#!/usr/bin/env python3
"""
This is a SeisFlows Test class which is used to test out the underlying
machinery before running an actual workflow. Contains simple functions used to
make sure that all parts of the package are working as expected.
"""
import os
import sys
import time
import logging
from glob import glob
from seisflows.tools import msg
from seisflows.config import (SeisFlowsPathsParameters, custom_import, ROOT_DIR,
                              CFGPATHS)

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]

system = sys.modules["seisflows_system"]
solver = sys.modules["seisflows_solver"]
optimize = sys.modules["seisflows_optimize"]
preprocess = sys.modules["seisflows_preprocess"]
postprocess = sys.modules["seisflows_postprocess"]


class Test(custom_import("workflow", "base")):
    """
    This is a template Base class
    """
    # Class-specific logger accessed using self.logger
    # When this logger is called, e.g., self.logger.info("text"), the logging
    # package will know exactly which module, class and function the log
    # statement has been sent from, extraordinarily helpful for debugging.
    logger = logging.getLogger(__name__).getChild(__qualname__)

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.

        :rtype: seisflows.config.SeisFlowsPathsParameters
        :return: Paths and parameters that define the given class

        """
        sf = SeisFlowsPathsParameters(super().required)

        sf.path("TEST_DATA", required=False,
                default=os.path.join(ROOT_DIR, "tests", "test_data",
                                     "work"),
                docstr="Example data for test system"
                )

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths. The validate function ensures that all
        required paths and parameters are accounted for, and that all
        optional paths and parameters are set to user-defined or default values.

        :type validate: bool
        :param validate: set required paths and parameters into sys.modules
        """
        if validate:
            self.required.validate()

    def main(self, return_flow=False):
        """
        This controls the main testing workflow
        """
        FLOW = [#self.test_system,
                # self.test_preprocess,
                self.test_solver,
                ]
        if return_flow:
            return FLOW

        for func in FLOW:
            func()

    def test_function(self, check_value):
        """
        A simple function that can be called by system.run()
        """
        print(f"Hello world, from taskid {system.taskid()}. "
              f"Check: {check_value}")

    def test_system(self):
        """
        Test the system by submitting a simple print statement using the
        run() and run(single) functions.

        Check that these functions perform as expected by passing in a random
        value and checking that this value gets logged back
        """
        check_value_1 = 1234.5
        system.run(classname="workflow", method="test_function",
                   check_value=check_value_1)

        time.sleep(3)  # wait a bit for system to catch up

        check_value_2 = 5432.1
        system.run(classname="workflow", method="test_function", single=True,
                   check_value=check_value_2)

        for fid, check in zip(
                sorted(glob(os.path.join(CFGPATHS.LOGDIR, "*.log"))),
                [check_value_1, check_value_2]
        ):
            with open(fid, "r") as f:
                line = f.readlines()[0]
                assert(float(line.strip().split(" ")[-1]) == check)

    def test_preprocess(self):
        """
        Test the exposed 'prepare_eval_grad()' preprocessing function
        """
        cwd = PATH.TEST_DATA
        taskid = 0
        filenames = ["AA.S0001.BXY.semd"]
        source_name = "001"
        preprocess.prepare_eval_grad(cwd=cwd, taskid=taskid,
                                     filenames=filenames,
                                     source_name=source_name
                                     )

    def test_solver(self):
        """
        Simply test that the solver binaries can be called, which is what the
        solver module is ultimately responsible for
        """
        assert(PATH.SPECFEM_BIN is not None and
               os.path.exists(PATH.SPECFEM_BIN)), (
            f"SPECFEM_BIN {PATH.SPECFEM_BIN} directory does not exist"
        )
        # SPECFEM2D won't have this executable
        try:
            solver.call_solver(
                executable=f"{PATH.SPECFEM_BIN}/xcombine_sem",
                output=os.path.join(PATH.TEST_DATA, "test_solver.log")
            )
        # We expect this to throw a system exit because
        except SystemExit:
            pass


    def test_postprocess(self):
        """

        """
        pass

