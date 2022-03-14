#!/usr/bin/env python
"""
This is the base class seisflows.workflow.Inversion

This is a main Seisflows class, it controls the main workflow.
"""
import os
import sys
import logging
import numpy as np
from glob import glob

from seisflows3.config import custom_import, CFGPATHS
from seisflows3.tools import msg, unix
from seisflows3.tools.wrappers import exists
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
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        super().__init__()

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

        # Define the Paths required by this module
        sf.path("MODEL_INIT", required=True,
                docstr="Initial model to be used for workflow")

        sf.path("MODEL_TRUE", required=False,
                docstr="Target model to be used for PAR.CASE == 'synthetic'")

        sf.path("DATA", required=False, default=None,
                docstr="path to data available to workflow")

        sf.path("FUNC", required=False,
                default=os.path.join(PATH.SCRATCH, CFGPATHS.SCRATCHDIR),
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
        msg.check(type(self))

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

    def main(self, return_flow=False):
        """
        !!! This function controls the main workflow !!!

        Carries out seismic inversion by running a series of functions in order

        :type return_flow: bool
        :param return_flow: for CLI tool, simply returns the flow function
            rather than running the workflow. Used for print statements etc.
        """

        # The workFLOW is a list of functions that can be called dynamically
        FLOW = [self.initialize,
                self.evaluate_gradient,
                self.write_gradient,
                self.compute_direction,
                self.line_search,
                self.finalize,
                self.clean
                ]
        if return_flow:
            return FLOW
        else:
            # FLOW is a constant, when we run, we flow 
            flow = FLOW

        optimize.iter = PAR.BEGIN

        # Allow workflow resume from a given mid-workflow location
        if PAR.RESUME_FROM:
            # Determine the index that corresponds to the resume function named
            try:
                resume_idx = [_.__name__ for _ in FLOW].index(PAR.RESUME_FROM)
                resume_fx = FLOW[resume_idx].__name__
            except ValueError:
                self.logger.info(f"{PAR.RESUME_FROM} does not correspond to any "
                                 f"workflow functions. Exiting...")
                sys.exit(-1)
            
            self.logger.info(
                    msg.mjr(f"RESUME ITERATION {optimize.iter} FROM FUNCTION: "
                            f"'{resume_fx}'")
                    )
            # Curtail the flow argument during resume_from
            flow = FLOW[resume_idx:]

        # First-time and one-time intialization of the workflow
        elif optimize.iter == 1:
            self.logger.info(msg.mjr("STARTING INVERSION WORKFLOW"))
            self.setup()

        # Run the workflow until from the current iteration until PAR.END
        while optimize.iter <= PAR.END:
            self.logger.info(msg.mnr(f"ITERATION {optimize.iter} / {PAR.END}"))

            # Execute the functions within the flow
            for func in flow:
                func()

                # Forcefully stop the workflow at STOP_AFTER (if requested)
                if PAR.STOP_AFTER and func.__name__ == PAR.STOP_AFTER:
                    self.logger.info(msg.mjr(f"STOP ITERATION {optimize.iter} "
                                             f"AT FUNCTION: '{PAR.STOP_AFTER}'")
                                     )
                    sys.exit(0)
            # Finish. Assuming completion of all arguments in flow() 
            self.logger.info(msg.mjr(f"FINISHED ITERATION {optimize.iter}"))
            optimize.iter += 1
            self.logger.info(msg.sub(f"SETTING ITERATION = {optimize.iter}"))
    
            # Reset flow incase 'resume_from' curtailed some of the arguments 
            flow = FLOW

    def setup(self):
        """
        Lays groundwork for inversion by running setup() functions for the 
        involved sub-modules, and generating synthetic true data if necessary, 
        and generating the pre-requisite database files. Should only be run once
        at the iteration 1
        """
        # Set up all the requisite modules from the master job
        self.logger.info(msg.mnr("PERFORMING MODULE SETUP"))
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
        self.logger.info(msg.mjr("INITIALIZING INVERSION"))
        self.evaluate_function(path=PATH.GRAD, suffix="new")

    def compute_direction(self):
        """
        Computes search direction
        """
        self.logger.info(msg.mnr("COMPUTING SEARCH DIRECTION"))
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
            self.logger.info(msg.mjr(f"CONDUCTING LINE SEARCH "
                                     f"({optimize.eval_str})")
                             )
            optimize.initialize_search()

        # Attempt a new trial step with the given step length
        optimize.line_search.step_count += 1
        self.logger.info(msg.mnr(f"TRIAL STEP COUNT: {optimize.eval_str}"))
        self.evaluate_function(path=PATH.FUNC, suffix="try")

        # Check the function evaluation against line search history
        status = optimize.update_search()

        # Proceed based on the outcome of the line search
        if status > 0:
            self.logger.info(msg.sub("trial step successful. finalizing..."))
            # Save outcome of line search to disk; reset step to 0 for next iter
            optimize.finalize_search()
            return
        elif status == 0:
            self.logger.info(msg.sub("retrying with new trial step"))
            # Recursively call this function to attempt another trial step
            self.line_search()
        elif status < 0:
            if optimize.retry_status():
                self.logger.info(msg.sub("line search failed. restarting..."))
                # Reset the line search machinery, do not reset step count (?)
                optimize.restart()
                self.line_search()
            else:
                self.logger.info(msg.sub("line search failed. aborting..."))
                sys.exit(-1)

    def evaluate_function(self, path, suffix):
        """
        Performs forward simulation, and evaluates the objective function

        :type path: str
        :param path: path in the scratch directory to use for I/O
        :type suffix: str
        :param suffix: suffix to use for I/O
        """
        self.logger.info(msg.sub("EVALUATE OBJECTIVE FUNCTION"))

        self.write_model(path=path, suffix=suffix)
        self.logger.debug(f"results saved to with suffix '{suffix}' to path: "
                          f"{path}")

        self.logger.info(f"evaluating objective function {PAR.NPROC} times")
        system.run("solver", "eval_func", path=path)

        self.logger.info("summing residuals with preprocess module")
        self.write_misfit(path=path, suffix=suffix)

    def evaluate_gradient(self, path=None):
        """
        Performs adjoint simulation to retrieve the gradient of the objective 
        """
        self.logger.info(msg.mnr("EVALUATING GRADIENT"))
        self.logger.info(f"evaluating gradient {PAR.NPROC} times")
        system.run("solver", "eval_grad", path=path or PATH.GRAD,
                   export_traces=PAR.SAVETRACES)

    def finalize(self):
        """
        Saves results from current model update iteration
        """
        self.logger.info(msg.mnr("FINALIZING WORKFLOW"))

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
        self.logger.info(msg.mnr("CLEANING WORKDIR FOR NEXT ITERATION"))

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

        self.logger.debug(f"saving model '{src}' to: {dst}")
        solver.save(solver.split(optimize.load(src)), dst)

    def write_gradient(self):
        """
        Writes gradient in format expected by non-linear optimization library.
        Calls the postprocess module, which will smooth/precondition gradient.
        """
        self.logger.info(msg.key("POSTPROCESSING KERNELS"))
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
        self.logger.debug(f"saving misfit {total_misfit:.3E} to: '{dst}'")
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

