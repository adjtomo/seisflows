#!/usr/bin/env python3
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

        # Define the Paths required by this module
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

    def main(self, return_flow=False):
        """
        This function controls the main SeisFlows3 workflow, and is submitted
        to system by the call `seisflows submit` or `seisflows resume`. It
        proceeds to evaluate a list of functions in order until a User defined
        stop criteria is met.

        :type return_flow: bool
        :param return_flow: for CLI tool, simply returns the flow function
            rather than running the workflow. Used for print statements etc.
        """
        # The workFLOW is a tuple of functions that can be called dynamic ally
        flow = (self.setup,
                self.initialize,
                self.evaluate_gradient,
                self.write_gradient,
                self.compute_direction,
                self.line_search,
                self.finalize,
                self.clean
                )
        if return_flow:
            return flow

        # Allow workflow resume from and stop after given flow functions
        start, stop = self.check_stop_resume_cond(flow)

        # Run the workflow until from the current iteration until PAR.END
        optimize.iter = PAR.BEGIN
        self.logger.info(msg.mjr("STARTING INVERSION WORKFLOW"))

        while True:
            self.logger.info(msg.mnr(f"ITERATION {optimize.iter} / {PAR.END}"))

            # Execute the functions within the flow
            for func in flow[start:stop]:
                func()

            # Finish. Assuming completion of all arguments in flow()
            self.logger.info(msg.mjr(f"FINISHED FLOW EXECUTION"))

            # Reset flow for subsequent iterations
            start, stop = None, None

            if optimize.iter == PAR.END:
                break

        self.logger.info(msg.mjr("FINISHED INVERSION WORKFLOW"))

    def setup(self):
        """
        Lays groundwork for inversion by running setup() functions for the 
        involved sub-modules, generating True model synthetic data if necessary,
        and generating the pre-requisite database files.

        .. note::
            This function should only be run one time, at the start of iter 1
        """
        # Iter check is done inside setup() so that we can include fx in FLOW
        if optimize.iter == 1:
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
            self.logger.info("trial step successful")
            # Save outcome of line search to disk; reset step to 0 for next iter
            optimize.finalize_search()
            return
        elif status == 0:
            self.logger.info("retrying with new trial step")
            # Recursively call this function to attempt another trial step
            self.line_search()
        elif status < 0:
            if optimize.retry_status():
                self.logger.info("line search failed. restarting line search")
                # Reset the line search machinery; set step count to 0
                optimize.restart()
                self.line_search()
            else:
                self.logger.info("line search failed. aborting inversion.")
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

        # Ensure that we are referencing the same tags as defined in OPTIMIZE
        model_tag = getattr(optimize, f"m_{suffix}")
        misfit_tag = getattr(optimize, f"f_{suffix}")

        self.write_model(path=path, tag=model_tag)

        self.logger.debug(f"evaluating objective function {PAR.NPROC} times")
        system.run("solver", "eval_func", path=path)

        self.write_misfit(path=path, tag=misfit_tag)

    def evaluate_gradient(self, path=None):
        """
        Performs adjoint simulation to retrieve the gradient of the objective 
        """
        self.logger.info(msg.mnr("EVALUATING GRADIENT"))
        self.logger.debug(f"evaluating gradient {PAR.NPROC} times")

        system.run("solver", "eval_grad", path=path or PATH.GRAD,
                   export_traces=PAR.SAVETRACES)

    def finalize(self):
        """
        Saves results from current model update iteration and increment the
        iteration number to set up for the next iteration. Finalization is
        expected to the be LAST function in workflow.main()'s  flow list.
        """
        self.logger.info(msg.mjr(f"FINALIZING ITERATION {optimize.iter}"))

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

        optimize.iter += 1
        self.logger.info(msg.sub(f"INCREMENT ITERATION TO {optimize.iter}"))

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

    def write_model(self, path, tag):
        """
        Writes model in format expected by solver

        :type path: str
        :param path: path to write the model to
        :type src: str
        :param src: name of the model to be saved, usually tagged as 'm' with
            a suffix depending on where in the inversion we are. e.g., 'm_try'.
            Expected that these tags are defined in OPTIMIZE module
        """
        src = tag
        dst = os.path.join(path, "model")
        self.logger.debug(f"saving model '{src}' to:\n{dst}")
        solver.save(solver.split(optimize.load(src)), dst)

    def write_gradient(self):
        """
        Writes gradient in format expected by non-linear optimization library.
        Calls the postprocess module, which will smooth/precondition gradient.
        """
        self.logger.info(msg.mnr("POSTPROCESSING KERNELS"))
        src = os.path.join(PATH.GRAD, "gradient")
        dst = f"g_new"

        postprocess.write_gradient(PATH.GRAD)
        parts = solver.load(src, suffix="_kernel")

        optimize.save(dst, solver.merge(parts))

    def write_misfit(self, path, tag):
        """
        Writes misfit in format expected by nonlinear optimization library.
        Collects all misfit values within the given residuals directory and sums
        them in a manner chosen by the preprocess class.

        :type path: str
        :param path: path to write the misfit to
        :type tag: str
        :param tag: name of the model to be saved, usually tagged as 'f' with
            a suffix depending on where in the inversion we are. e.g., 'f_try'.
            Expected that these tags are defined in OPTIMIZE module
        """
        self.logger.info("summing residuals with preprocess module")
        src = glob(os.path.join(path, "residuals", "*"))
        dst = tag
        total_misfit = preprocess.sum_residuals(src)

        self.logger.debug(f"saving misfit {total_misfit:.3E} to tag '{dst}'")
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
            src = os.path.join(PATH.OPTIMIZE, optimize.g_old)
            unix.cp(src, dst + ".npy")

        self.logger.debug(f"saving gradient to path:\n{dst}")

    def save_model(self):
        """
        Save the model vector. Allows saving numpy array or standard
        Fortran .bin files

        Saving as a vector saves on file count, but requires numpy and seisflows
        functions to read
        """
        src = optimize.m_new
        dst = os.path.join(PATH.OUTPUT, f"model_{optimize.iter:04d}")

        self.logger.debug(f"saving model '{src}' to path:\n{dst}")

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

        self.logger.debug(f"saving kernels to path:\n{dst}")

        unix.mv(src, dst)

    def save_traces(self):
        """
        Save the waveform traces to disk.

        !!! This doesn't work? Traces are not saved to PATH.GRAD so src does
        !!! not exist
        """
        src = os.path.join(PATH.GRAD, "traces")
        dst = os.path.join(PATH.OUTPUT, f"traces_{optimize.iter:04d}")

        self.logger.debug(f"saving traces to path:\n{dst}")

        unix.mv(src, dst)

    def save_residuals(self):
        """
        Save the residuals to disk
        """
        src = os.path.join(PATH.GRAD, "residuals")
        dst = os.path.join(PATH.OUTPUT, f"residuals_{optimize.iter:04d}")

        self.logger.debug(f"saving residuals to path:\n{dst}")

        unix.mv(src, dst)

