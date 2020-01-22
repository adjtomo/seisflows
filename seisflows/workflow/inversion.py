#!/usr/bin/env python
"""
This is the base class seisflows.workflow.Inversion

This is a main Seisflows class, it controls the main workflow.
"""
import os
import sys
import time
from glob import glob

import numpy as np
from seisflows.tools import unix
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
postprocess = sys.modules["seisflows_postprocess"]


class Inversion(Base):
    """
    Waveform inversion base class

    Peforms iterative nonlinear inversion and provides a base class on top
    of which specialized strategies can be implemented.

    To allow customization, the inversion workflow is divided into generic
    methods such as "initialize", "finalize", "evaluate_function",
    "evaluate_gradient", which can be easily overloaded.

    Calls to forward and adjoint solvers are abstracted through the "solver"
    interface so that various forward modeling packages can be used
    interchangeably.

    Commands for running in serial or parallel on a workstation or cluster
    are abstracted through the "system" interface.
    """
    def check(self):
        """
        Checks parameters and paths
        """
        # Starting and stopping iterations
        if "BEGIN" not in PAR:
            raise ParameterError(PAR, "BEGIN")

        if "END" not in PAR:
            raise ParameterError(PAR, "END")

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

        # Input paths
        if "DATA" not in PATH:
            setattr(PATH, "DATA", None)

        if "MODEL_INIT" not in PATH:
            raise ParameterError(PATH, "MODEL_INIT")

        # Output paths
        if "OUTPUT" not in PATH:
            raise ParameterError(PATH, "OUTPUT")

        # Permanet disk exports
        if "SAVEMODEL" not in PAR:
            setattr(PAR, "SAVEMODEL", True)

        if "SAVEGRADIENT" not in PAR:
            setattr(PAR, "SAVEGRADIENT", False)

        if "SAVEKERNELS" not in PAR:
            setattr(PAR, "SAVEKERNELS", False)

        if "SAVEAS" not in PAR:
            setattr(PAR, "SAVEAS", "binary")

        if "SAVETRACES" not in PAR:
            setattr(PAR, "SAVETRACES", False)

        if "SAVERESIDUALS" not in PAR:
            setattr(PAR, "SAVERESIDUALS", False)

        # Print statement outputs
        if "VERBOSE" not in PAR:
            setattr(PAR, "VERBOSE", True)

        # Parameter assertions
        assert 1 <= PAR.BEGIN <= PAR.END

        # Path assertions
        if not exists(PATH.DATA):
            assert "MODEL_TRUE" in PATH
            assert exists(PATH.MODEL_TRUE)

        # Check that there is a given starting model
        if not exists(PATH.MODEL_INIT):
            raise Exception("MODEL_INIT does not exist")

    def stopwatch(self, method):
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
        !!! This function controls the main workflow!!!

        Carries out seismic inversion
        """
        self.stopwatch("set")
        print(f"BEGINNING WORKFLOW AT {self.stopwatch()}")

        # One-time intialization of the workflow
        optimize.iter = PAR.BEGIN
        self.setup()
        print(f"{optimize.iter} <= {PAR.END}")

        # Run the workflow until PAR.END
        while optimize.iter <= PAR.END:
            print(f"ITERATION {optimize.iter}")
            self.initialize()
            self.evaluate_gradient()
            self.compute_direction()
            self.line_search()
            self.finalize()
            self.clean()
            print(f"finished iteration {optimize.iter} at {self.stopwatch()}\n")
            self.stopwatch("time")
            optimize.iter += 1

    def setup(self):
        """
        Lays groundwork for inversion
        """
        # Set up all the requisite modules
        if optimize.iter == 1:
            print("SETUP\n\tPerforming module setup")
            preprocess.setup()
            postprocess.setup()
            optimize.setup()

        if optimize.iter == 1 or PATH.LOCAL:
            if PATH.DATA:
                print("Copying data")
            else:
                print("Generating data")

            print("Running solver")
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

    def compute_direction(self):
        """
        Computes search direction
        """
        print("COMPUTE SEARCH DIRECTION")
        optimize.compute_direction()

    def line_search(self):
        """ Conducts line search in given search direction

          Status codes
              status > 0  : finished
              status == 0 : not finished
              status < 0  : failed
        """
        print("LINE SEARCH\n\tinitializing line search")
        optimize.initialize_search()

        while True:
            print(f"TRIAL STEP: {optimize.line_search.step_count + 1}")
            self.evaluate_function()
            status = optimize.update_search()
            # Determine the outcome of the line search
            if status > 0:
                print("\tTrial step successful")
                optimize.finalize_search()
                break
            elif status == 0:
                print("\tRetrying with new trial step")
                continue
            elif status < 0:
                if optimize.retry_status():
                    print("\tLine search failed\n\n Retrying...")
                    optimize.restart()
                    self.line_search()
                    break
                else:
                    print("\tLine search failed\n\n Aborting...")
                    sys.exit(-1)

    def evaluate_function(self):
        """
        Performs forward simulation to evaluate objective function
        """
        self.write_model(path=PATH.FUNC, suffix="try")

        system.run("solver", "eval_func", path=PATH.FUNC)

        self.write_misfit(path=PATH.FUNC, suffix="try")

    def evaluate_gradient(self):
        """
        Performs adjoint simulation to evaluate gradient of the misfit function
        """
        print("EVALUATE GRADIENT\n\tRunning adjoint simulation")

        self.stopwatch("set")
        system.run("solver", "eval_grad", path=PATH.GRAD,
                   export_traces=PAR.SAVETRACES)
        self.stopwatch("time")
        self.write_gradient(path=PATH.GRAD, suffix="new")

    def finalize(self):
        """
        Saves results from current model update iteration
        """
        self.checkpoint()

        # Save files from scratch before discarding
        if PAR.SAVEMODEL:
            self.save_model()

        if PAR.SAVEGRADIENT:
            self.save_gradient()

        if PAR.SAVEKERNELS:
            self.save_kernels()

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVERESIDUALS:
            self.save_residuals()

    def clean(self):
        """
        Cleans directories in which function and gradient evaluations were
        carried out
        """
        print("CLEAN")

        unix.rm(PATH.GRAD)
        unix.rm(PATH.FUNC)
        unix.mkdir(PATH.GRAD)
        unix.mkdir(PATH.FUNC)

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

    def write_gradient(self, path="", suffix=""):
        """
        Writes gradient in format expected by nonlinear optimization library

        :type path: str
        :param path: path to write the model to
        :type suffix: str
        :param suffix: suffix to add to the model
        """
        print("POSTPROCESSING")
        src = os.path.join(path, "gradient")
        dst = f"g_{suffix}"

        self.stopwatch("set")
        postprocess.write_gradient(path)
        parts = solver.load(src, suffix="_kernel")
        optimize.save(dst, solver.merge(parts))
        self.stopwatch("time")

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

    def create_vtk_file(self, src, dst):
        src = os.path.join(PATH.GRAD)
        raise NotImplementedError

    def save_gradient(self):
        """
        Save the gradient vector. Allows saving numpy array or standard
        Fortran .bin files

        Saving as a vector saves on file count, but requires numpy and seisflows
        functions to read
        """
        dst = os.path.join(PATH.OUTPUT, f"gradient_{optimize.iter:04d}")

        if PAR.SAVEAS in ["binary", "both"]:
            src = os.path.join(PATH.GRAD, "gradient")
            unix.mv(src, dst)
        if PAR.SAVEAS in ["vector", "both"]:
            src = os.path.join(PATH.OPTIMIZE, "g_old")
            unix.cp(src, dst + ".npy")

    def save_model(self):
        """
        Save the model vector. Allows saving numpy array or standard
        Fortran .bin files

        Saving as a vector saves on file count, but requires numpy and seisflows
        functions to read
        """
        src = "m_new"
        dst = os.path.join(PATH.OUTPUT, "model_{optimize.iter:04d}")
        if PAR.SAVEAS in ["binary", "both"]:
            solver.save(solver.split(optimize.load(src)), dst)
        if PAR.SAVEAS in ["vector", "both"]:
            np.save(file=dst, arr=optimize.load(src))

    def save_kernels(self):
        """
        Save the kernel vector as a Fortran binary file
        """
        src = os.path.join(PATH.GRAD, "kernels")
        dst = os.path.join(PATH.OUTPUT, "kernels_{optimize.iter:04d}")
        unix.mv(src, dst)

    def save_traces(self):
        """
        Save the traces
        """
        src = os.path.join(PATH.GRAD, "traces")
        dst = os.path.join(PATH.OUTPUT, "traces_{optimize.iter:04d}")
        unix.mv(src, dst)

    def save_residuals(self):
        """
        Save the residuals
        """
        src = os.path.join(PATH.GRAD, "residuals")
        dst = os.path.join(PATH.OUTPUT, "residuals_{optimize.iter:04d}")
        unix.mv(src, dst)

