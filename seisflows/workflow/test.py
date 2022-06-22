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
import subprocess
import numpy as np

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
                # self.test_solver,
                self.test_optimize
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
        # Run a very simple test function using system.run()
        check_value_1 = 1234.5
        system.run(classname="workflow", method="test_function",
                   check_value=check_value_1)

        time.sleep(3)  # wait a bit for system to catch up

        check_value_2 = 5432.1
        system.run(classname="workflow", method="test_function", single=True,
                   check_value=check_value_2)

        # Check the output log files to match the check values
        for fid, check in zip(
                sorted(glob(os.path.join(CFGPATHS.LOGDIR, "*.log"))),
                [check_value_1, check_value_2]
        ):
            with open(fid, "r") as f:
                line = f.readlines()[0]
                assert(float(line.strip().split(" ")[-1]) == check)

        # Check that MPI Exec works
        assert("MPIEXEC" in PAR), f"MPIEXEC is not defined for this system"
        stdout = subprocess.run(PAR.MPIEXEC, shell=True, check=True,
                                stdout=subprocess.PIPE)

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
        try:
            solver.call_solver(
                executable=f"{PATH.SPECFEM_BIN}/xcombine_sem",
                output=os.path.join(PATH.TEST_DATA, "test_solver.log")
            )
        # We expect this to throw a system exit because we are not running with
        # MPI
        except SystemExit:
            pass

    def test_optimize(self):
        """
        Test optimization module with a simple Rosenbrock function
        """
        m_new, m_true, objective_function, gradient = rosenbrock()
        optimize.setup(m_new=m_new)

        def evaluate_function():
            """
            Evalaute the misfit function of a given model
            """
            self.logger.info("evaluating objective function")
            m_try = optimize.load("m_try")
            f_try = objective_function(m_try)
            optimize.save("f_try", f_try)

        def evaluate_gradient():
            """
            Evaluate the gradient of a given model
            """
            self.logger.info("evaluating gradient")
            m_new = optimize.load("m_new")
            f_new = objective_function(m_new)
            g_new = gradient(m_new)
            optimize.save("f_new", f_new)
            optimize.save("g_new", g_new)

        def line_search():
            """
            Run a line search until a suitable model has been found
            Note this is almost the same as workflow.inversion.line_search
            """
            optimize.initialize_search()
            while True:
                evaluate_function()
                optimize.line_search.step_count += 1
                status = optimize.update_search()
                if status == 1:
                    self.logger.info("finalizing line search")
                    optimize.finalize_search()
                    return
                elif status == 0:
                    self.logger.info("continuing line search")
                    continue
                elif status == -1:
                    if optimize.retry_status():
                        self.logger.info("restarting line search")
                        optimize.restart()
                        # Recursively run the line search after restart
                        line_search()
                    else:
                        sys.exit(-1)

        def finalize(thresh=5e-3):
            """
            Finish off one iteration, check the distance between old and new
            vectors to see if model stops changing
            """
            m_new = optimize.load("m_new")
            m_diff = np.linalg.norm(m_new - m_true) / np.linalg.norm(m_new)
            if m_diff < thresh:
                self.logger.info(f"successful inversion after {optimize.iter} "
                                 f"iterations")
                sys.exit(0)
            else:
                self.logger.info(f"model difference: {m_diff:.2E}")
                return

        self.logger.info("testing optimization library with Rosenbrock problem")
        for iteration in range(1, 200):
            self.logger.info(f"iteration {iteration}")
            evaluate_gradient()
            optimize.compute_direction()
            line_search()
            optimize.iter += 1
            finalize()


def rosenbrock():
    """
    Rosenbrock test problem for optimization library testing

    https://en.wikipedia.org/wiki/Rosenbrock_function
    """
    model_init = np.array([-1.2, 1])  # This is the guess for the global min
    model_true = np.array([1, 1])  # This is the actual minimum

    def objective_function(x):
        """
        Rosenbrock objective function which is defined mathematically as:

        f(x,y) = (a-x)^2 + b(y-x^2)^2

        where the global minimum is at (x,y) == (a, a^2)
        and typical constant values are: a==1, b==100
        """
        return np.array([((1 - x[0]) ** 2 + 100 * (-x[0] ** 2 + x[1]) ** 2)])

    def gradient(x):
        """
        Gradient of the objective function for Rosenbrock test
        """
        return np.array([-2*(1-x[0]) - 400*x[0]*(-x[0]**2+x[1]),
                         200*(- x[0]**2+x[1])])

    return model_init, model_true, objective_function, gradient


def rosenbrock_n(n=1E5):
    """
    N dimensional Rosenbrock test problem for optimization library testing

    https://en.wikipedia.org/wiki/Rosenbrock_function
    """
    model_init = 0.1 * np.ones(int(n))  # This is a guess for the global min
    model_true = np.ones(int(n))

    def objective_function(x):
        """
        Rosenbrock objective function
        """
        return sum(100 * (x[:-1]**2. - x[1:])**2. + (x[:-1] - 1.)**2)

    def gradient(x):
        """
        Gradient of the objective function for Rosenbrock test
        """
        g = np.zeros(int(n))
        g[1:-1] = -200 * (x[:-2] ** 2. - x[1:-1]) + \
                  400. * x[1:-1] * (x[1:-1] ** 2. - x[2:]) + \
                  2. * (x[1:-1] - 1.)

        g[0] = 400. * x[0] * (x[0] ** 2. - x[1]) + \
               2. * (x[0] - 1)

        g[-1] = -200. * (x[-2] ** 2. - x[-1])

        return g

    return model_init, model_true, objective_function, gradient

