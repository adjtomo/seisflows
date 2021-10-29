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
from seisflows3.config import custom_import
from seisflows3.tools import unix
from seisflows3.tools.tools import exists
from seisflows3.config import save, SeisFlowsPathsParameters

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
    interface so that various forward modeling packages canf be used
    interchangeably.

    Commands for running in serial or parallel on a workstation or cluster
    are abstracted through the "system" interface.
    """
    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        # Define the Parameters required by this module
        sf.par("BEGIN", required=True, par_type=int,
               docstr="First iteration of workflow, 1 <= BEGIN <= inf")

        sf.par("END", required=True, par_type=int,
               docstr="Last iteration of workflow, BEGIN <= END <= inf")

        sf.par("RESUME_FROM", required=False, par_type=str,
               docstr="Name of task to resume inversion from")

        sf.par("STOP_AFTER", required=False, par_type=str,
               docstr="Name of task to stop inversion after finishing")

        sf.par("CASE", required=True, par_type=str,
               docstr="Type of inversion, available: "
                      "['data': real data inversion, "
                      "'synthetic': synthetic-synthetic inversion]")

        sf.par("SAVEMODEL", required=False, default=True, par_type=bool,
               docstr="Save final model files after each iteration")

        sf.par("SAVEGRADIENT", required=False, default=True, par_type=bool,
               docstr="Save gradient files after each iteration")

        sf.par("SAVEKERNELS", required=False, default=False, par_type=bool,
               docstr="Save event kernel files after each iteration")

        sf.par("SAVETRACES", required=False, default=False, par_type=bool,
               docstr="Save waveform traces after each iteration")

        sf.par("SAVERESIDUALS", required=False, default=False, par_type=bool,
               docstr="Save waveform residuals after each iteration")

        sf.par("SAVEAS", required=False, default="binary", par_type=str,
               docstr="Format to save models, gradients, kernels. "
                      "Available: "
                      "['binary': save files in native SPECFEM .bin format, "
                      "'vector': save files as NumPy .npy files, "
                      "'both': save as both binary and vectors]")

        sf.par("VERBOSE", required=False, default=True, par_type=bool,
               docstr="Provide detailed statements to the output logs")

        # Define the Paths required by this module
        sf.path("MODEL_INIT", required=True,
                docstr="Initial model to be used for workflow")

        sf.path("MODEL_TRUE", required=False,
                docstr="Target model to be used for PAR.CASE == 'synthetic'")

        sf.path("DATA", required=False, default=None,
                docstr="path to data available to workflow")

        sf.path("FUNC", required=False,
                default=os.path.join(PATH.SCRATCH, "evalfunc"),
                docstr="scratch path to store data related to function "
                       "evaluations")

        sf.path("GRAD", required=False,
                default=os.path.join(PATH.SCRATCH, "evalgrad"),
                docstr="scratch path to store data related to gradient "
                       "evaluations")

        sf.path("HESS", required=False,
               default=os.path.join(PATH.SCRATCH, "evalhess"),
               docstr="scratch path to store data related to Hessian "
                      "evaluations")

        sf.path("OPTIMIZE", required=False,
                default=os.path.join(PATH.SCRATCH, "optimize"),
                docstr="scratch path to store data related to nonlinear "
                       "optimization")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=False)
        if validate:
            self.required.validate()

        for required_path in ["SCRATCH", "OUTPUT", "LOCAL"]:
            assert(required_path in PATH), \
                f"Inversion requires path {required_path}"

        assert(1 <= PAR.BEGIN <= PAR.END), \
            f"Incorrect BEGIN or END parameter: 1 <= {PAR.BEGIN} <= {PAR.END}"

        if PAR.CASE.upper() == "SYNTHETIC":
            assert exists(PATH.MODEL_TRUE), \
                "CASE == SYNTHETIC requires PATH.MODEL_TRUE"

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
        # Run solver.setup() in parallel
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
        # Calculate the initial step length based on optimization algorithm
        if optimize.line_search.step_count == 0:
            print("LINE SEARCH")
            optimize.initialize_search()

        # Attempt a new trial step with the given step length
        optimize.line_search.step_count += 1
        print(f"TRIAL STEP: {optimize.line_search.step_count}")
        self.evaluate_function(path=PATH.FUNC, suffix="try")

        # Check the function evaluation against line search history
        status = optimize.update_search()

        # Proceed based on the outcome of the line search
        if status > 0:
            print("\tTrial step successful")
            # Save outcome of line search to disk; reset step to 0 for next iter
            optimize.finalize_search()
            return
        elif status == 0:
            print("\tRetrying with new trial step")
            # Recursively call this function to attempt another trial step
            self.line_search()
        elif status < 0:
            if optimize.retry_status():
                print("\tLine search failed\n\n Restarting optimization...")
                # Reset the line search machinery, do not reset step count (?)
                optimize.restart()
                self.line_search()
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
        print("EVALUATE FUNCTION\n\tRunning forward simulation")
        self.write_model(path=path, suffix=suffix)
        system.run("solver", "eval_func", path=path)
        self.write_misfit(path=path, suffix=suffix)

    def evaluate_gradient(self, path=None):
        """
        Performs adjoint simulation to retrieve the gradient of the objective 
        """
        print("EVALUATE GRADIENT\n\tRunning adjoint simulation")
        system.run("solver", "eval_grad", path=path or PATH.GRAD,
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

