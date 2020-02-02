#!/usr/bin/env python
"""
This is the base class seisflows.workflow.forward

This is a main Seisflows class, it controls the main workflow.
"""
import os
import sys
import time
from glob import glob

from seisflows.tools.tools import exists
from seisflows.config import save
from seisflows.workflow.base import Base
from seisflows.tools.err import ParameterError

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]

system = sys.modules["seisflows_system"]
solver = sys.modules["seisflows_solver"]
optimize = sys.modules["seisflows_optimize"]
preprocess = sys.modules["seisflows_preprocess"]


class Forward(Base):
    """
    Waveform forward solver base class, performs a forward simulation for
    NTASK events. This is useful for model appraisal, misfit quanitification
    or for using the machinery of Seisflows to batch run forward simulations
    """

    def check(self):
        """
        Checks parameters and paths
        """
        # Scratch paths
        if "SCRATCH" not in PATH:
            raise ParameterError(PATH, "SCRATCH")

        if "LOCAL" not in PATH:
            setattr(PATH, "LOCAL", None)

        if "FUNC" not in PATH:
            setattr(PATH, "FUNC", os.path.join(PATH.SCRATCH, "evalfunc"))

        if "GRAD" not in PATH:
            setattr(PATH, "GRAD", os.path.join(PATH.SCRATCH, "evalgrad"))

        if "HESS" not in PATH:
            setattr(PATH, "HESS", os.path.join(PATH.SCRATCH, "evalhess"))

        if "OPTIMIZE" not in PATH:
            setattr(PATH, "OPTIMIZE", os.path.join(PATH.SCRATCH, "optimize"))

        if "FUNC" not in PATH:
            setattr(PATH, "FUNC", os.path.join(PATH.SCRATCH, "evalfunc"))

        # Input paths
        if "DATA" not in PATH:
            setattr(PATH, "DATA", None)

        if "MODEL_INIT" not in PATH:
            raise ParameterError(PATH, "MODEL_INIT")

        # Output paths
        if "OUTPUT" not in PATH:
            raise ParameterError(PATH, "OUTPUT")

        # Print statement outputs
        if "VERBOSE" not in PAR:
            setattr(PAR, "VERBOSE", True)

        # Parameter assertions
        assert 1 <= PAR.BEGIN <= PAR.END

        # Path assertions
        if PATH.DATA in PAR and not exists(PATH.DATA):
            assert "MODEL_TRUE" in PATH
            assert exists(PATH.MODEL_TRUE)

        # Check that there is a given starting model
        if not exists(PATH.MODEL_INIT):
            raise Exception("MODEL_INIT does not exist")

    def stopwatch(self, method=None):
        """
        Timestamps for print statements. `time` package wrapper

        :type method: str
        :param method: interact with the stopwatch:
            "set" to start the timer, "time" to return time since `set` time
            or None to return the current time
        :rtype: time.time
        :return: a time object
        """
        if method == "set" or not self._starttime:
            self._starttime = time.time()
        elif method == "time":
            elapsed = (time.time() - self._starttime) / 60.
            print(f"{elapsed:.2f}m elapsed")
        else:
            return time.asctime()

    def main(self):
        """
        !!! This function controls the main workflow !!!

        Carries out seismic inversion
        """
        self.stopwatch("set")
        print(f"BEGINNING WORKFLOW AT {self.stopwatch()}")
        self.setup()
        self.initialize()
        print(f"FINISHED AT {self.stopwatch()}\n")
        self.stopwatch("time")

    def setup(self):
        """
        Lays groundwork for inversion
        """
        preprocess.setup()
        print("Running solver", end="... ")
        self.stopwatch("set")
        system.run("solver", "setup")
        self.stopwatch("time")

    def initialize(self):
        """
        Generates synthetics via a forward simulation, calculates misfits
        for the forward simulation. Writes misfit for use in optimization.
        """
        self.write_model(path=PATH.GRAD, suffix="new")

        print("Generating synthetics")
        system.run("solver", "eval_func", path=PATH.GRAD)

        self.write_misfit(suffix="new")

    def checkpoint(self):
        """
        Writes information to disk so workflow can be resumed following a break
        """
        save()

    def write_model(self, path="", suffix=""):
        """
        Writes model in format expected by solver

        :type path: str
        :param path: path to write the model to
        :type suffix: str
        :param suffix: suffix to add to the model
        """
        src = f"m_{suffix}"
        dst = os.path.join(path, "model")

        solver.save(solver.split(optimize.load(src)), dst)

    def write_misfit(self, path="", suffix=""):
        """
        Writes misfit in format expected by nonlinear optimization library

        :type path: str
        :param path: path to write the misfit to
        :type suffix: str
        :param suffix: suffix to add to the misfit
        """
        src = glob(os.path.join(path, "residuals", "*"))
        dst = "f_{suffix}"

        total_misfit = preprocess.sum_residuals(src)
        optimize.savetxt(dst, total_misfit)
