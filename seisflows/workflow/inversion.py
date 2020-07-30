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
from seisflows.config import custom_import
from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import save
from seisflows.tools.err import ParameterError

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]

system = sys.modules["seisflows_system"]
solver = sys.modules["seisflows_solver"]
optimize = sys.modules["seisflows_optimize"]
preprocess = sys.modules["seisflows_preprocess"]
postprocess = sys.modules["seisflows_postprocess"]


class Inversion(custom_import("workflow", "base")):
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

        # Outputs to disk 
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
        if PATH.DATA in PAR and not exists(PATH.DATA):
            assert "MODEL_TRUE" in PATH, "MODEL_TRUE must be in PATH"
            assert exists(PATH.MODEL_TRUE), "MODEL_TRUE does not exist"

        # Check that there is a given starting model
        if not exists(PATH.MODEL_INIT):
            raise Exception("MODEL_INIT does not exist")

        # Check if this is a synthetic-synthetic or data-synthetic inversion
        if "CASE" not in PAR:
            raise ParameterError(PAR, "CASE")
        elif PAR.CASE.upper() == "SYNTHETIC" and not exists(PATH.MODEL_TRUE):
            raise Exception("CASE == SYNTHETIC requires PATH.MODEL_TRUE")

        if "RESUME_FROM" not in PAR:
            setattr(PAR, "RESUME_FROM", None)

        if "STOP_AFTER" not in PAR:
            setattr(PAR, "STOP_AFTER", None)

    def main(self):
        """
        !!! This function controls the main workflow !!!

        Carries out seismic inversion by running a series of functions in order
        """
        # The workflow is a list of functions that can be called dynamically
        flow = [self.initialize,
                self.evaluate_gradient,
                self.write_gradient,
                self.compute_direction,
                self.line_search,
                self.finalize,
                self.clean
                ]

        print(f"BEGINNING WORKFLOW AT {time.asctime()}")
        optimize.iter = PAR.BEGIN
        print(f"{optimize.iter} <= {PAR.END}")

        # Allow workflow resume from a given mid-workflow location
        if PAR.RESUME_FROM:
            self.resume_from(flow)
        elif optimize.iter == 1:
            # First-time intialization of the workflow
            self.setup()

        # Run the workflow until from the current iteration until PAR.END
        while optimize.iter <= PAR.END:
            print(f"ITERATION {optimize.iter}")
            for func in flow:
                func()
                # Stop the workflow at STOP_AFTER if requested
                if PAR.STOP_AFTER and func.__name__ == PAR.STOP_AFTER:
                    print(f"STOP ITERATION {optimize.iter} AT {PAR.STOP_AFTER}")
                    break
            print(f"FINISHED ITERATION {optimize.iter} AT {time.asctime()}\n")
            optimize.iter += 1

    def resume_from(self, flow):
        """
        Resume the workflow from a given function, proceed in the same fashion 
        as main until the end of the current iteration.

        :type flow: list of functions
        :param flow: the order of functions defined in main(), which should be
            parsed to determine where workflow should be resumed from
        """ 
        # Determine the index that corresponds to the resume function named
        try:
            resume_idx = [_.__name__ for _ in flow].index(PAR.RESUME_FROM)
        except ValueError:
            print(f"{PAR.RESUME_FROM} does not correspond to any workflow "
                  "functions. Exiting...")
            sys.exit(-1)

        print(f"RESUME ITERATION {optimize.iter} (from function "
              f"{flow[resume_idx].__name__})")
        
        for func in flow[resume_idx:]:
            func()

        print(f"FINISHED ITERATION {optimize.iter} AT {time.asctime()}\n")
        optimize.iter += 1

    def setup(self):
        """
        Lays groundwork for inversion by running setup() functions for the 
        involved sub-modules, and generating synthetic true data if necessary, 
        and generating the pre-requisite database files. Should only be run once
        at the iteration 1
        """
        # Set up all the requisite modules
        print("SETUP")
        preprocess.setup()
        postprocess.setup()
        optimize.setup()
        system.run("solver", "setup")

    def initialize(self):
        """
        Generates synthetics via a forward simulation, calculates misfits
        for the forward simulation. Writes misfit for use in optimization.
        """
        print("INITIALIZE")
        self.evaluate_function(path=PATH.GRAD, suffix="new")

    def compute_direction(self):
        """
        Computes search direction
        """
        print("COMPUTE SEARCH DIRECTION")
        optimize.compute_direction()

    def line_search(self):
        """
        Conducts line search in given search direction
        Status codes:
            status > 0  : finished
            status == 0 : not finished
            status < 0  : failed
        """
        print("LINE SEARCH")
        optimize.initialize_search()

        while True:
            print(f"TRIAL STEP: {optimize.line_search.step_count + 1}")
            self.evaluate_function(path=PATH.FUNC, suffix="try")
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
                    print("\tLine search failed\n\n Restarting optimization...")
                    optimize.restart()
                    self.line_search()
                    break
                else:
                    print("\tLine search failed\n\n Aborting...")
                    sys.exit(-1)

    def evaluate_function(self, path, suffix):
        """
        Performs forward simulation, and evaluates the objective function

        :type path: str
        :param path: path in the scratch directory to use for I/O
        :type suffix: str
        :param suffix: suffix to use for I/O
        """
        print("EVALUATION FUNCTION\n\tRunning forward simulation")
        self.write_model(path=path, suffix=suffix)
        system.run("solver", "eval_func", path=path)
        self.write_misfit(path=path, suffix=suffix)

    def evaluate_gradient(self):
        """
        Performs adjoint simulation to retrieve the gradient of the objective 
        """
        print("EVALUATE GRADIENT\n\tRunning adjoint simulation")
        system.run("solver", "eval_grad", path=PATH.GRAD,
                   export_traces=PAR.SAVETRACES)

    def finalize(self):
        """
        Saves results from current model update iteration
        """
        self.checkpoint()
        preprocess.finalize()

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

    def write_model(self, path, suffix):
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

    def write_gradient(self):
        """
        Writes gradient in format expected by non-linear optimization library.
        Calls the postprocess module, which will smooth/precondition gradient.
        """
        print("POSTPROCESSING")
        src = os.path.join(PATH.GRAD, "gradient")
        dst = f"g_new"

        postprocess.write_gradient(PATH.GRAD)
        parts = solver.load(src, suffix="_kernel")
        optimize.save(dst, solver.merge(parts))

    def write_misfit(self, path, suffix):
        """
        Writes misfit in format expected by nonlinear optimization library.
        Collects all misfit values within the given residuals directory and sums
        them in a manner chosen by the preprocess class.

        :type path: str
        :param path: path to write the misfit to
        :type suffix: str
        :param suffix: suffix to add to the misfit
        """
        src = glob(os.path.join(path, "residuals", "*"))
        dst = f"f_{suffix}"

        total_misfit = preprocess.sum_residuals(src)
        optimize.savetxt(dst, total_misfit)

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
        dst = os.path.join(PATH.OUTPUT, f"model_{optimize.iter:04d}")
        if PAR.SAVEAS in ["binary", "both"]:
            solver.save(solver.split(optimize.load(src)), dst)
        if PAR.SAVEAS in ["vector", "both"]:
            np.save(file=dst, arr=optimize.load(src))

    def save_kernels(self):
        """
        Save the kernel vector as a Fortran binary file on disk
        """
        src = os.path.join(PATH.GRAD, "kernels")
        dst = os.path.join(PATH.OUTPUT, f"kernels_{optimize.iter:04d}")
        unix.mv(src, dst)

    def save_traces(self):
        """
        Save the waveform traces to disk
        """
        src = os.path.join(PATH.GRAD, "traces")
        dst = os.path.join(PATH.OUTPUT, f"traces_{optimize.iter:04d}")
        unix.mv(src, dst)

    def save_residuals(self):
        """
        Save the residuals to disk
        """
        src = os.path.join(PATH.GRAD, "residuals")
        dst = os.path.join(PATH.OUTPUT, f"residuals_{optimize.iter:04d}")
        unix.mv(src, dst)

