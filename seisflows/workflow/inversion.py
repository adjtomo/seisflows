#!/usr/bin/env python3
"""
A seismic inversion (a.k.a full waveform inversion, adjoint tomography, full
waveform tomography) perturbs seismic velocity models by minimizing objective
functions defining differences between observed and synthetic waveforms.

This seismic inversion workflow performs a linear set of tasks involving:

1) Generating synthetic seismograms using an external numerical solver
2) Calculating time-dependent misfit (adjoint sources) between data
    (or other synthetics) and synthetics
3) Using adjoint sources to generate misfit kernels defining volumetric
    perturbations sensitive to data-synthetic misfit
4) Smoothing and summing misfit kernels into a single gradient
5) Perturbing the starting model with the gradient to reduce misfit defined by
    the objective function during a line search

The Inversion workflow runs the above tasks in a loop (iterations) while
exporting updated models, kernels and/or gradients to disk.
"""
import os
import sys
import numpy as np

from seisflows.workflow.migration import Migration
from seisflows.tools import msg, unix


class Inversion(Migration):
    """
    Waveform inversion base class, built on top of the Migration child class,
    which in-turn is built on top of the Forward child class.

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
    def __init__(self):
        """
        Init is used to instantiate global parameters defined by the input
        parameter file.
        """
        super().__init__()

        self.required.par(
            "BEGIN", required=False, default=1, par_type=int,
            docstr="First iteration of an inversion workflow, 1 <= BEGIN <= inf"
        )

        self.required.par(
            "END", required=False, default=1, par_type=int,
            docstr="Last iteration of the inverison workflow,"
                   "BEGIN <= END <= inf"
        )
        self.required.par(
            "RESUME_FROM", required=False, par_type=str,
            docstr="Name of flow task to resume workflow from. Useful for "
                   "restarting failed workflows or re-trying sections of "
                   "workflows with new parameters. To determine available "
                   "options for your given workflow: > seisflows print flow"
        )
        self.required.par(
            "STOP_AFTER", required=False, par_type=str,
            docstr="Name of flow task to stop workflow after. Useful for "
                   "stopping mid-workflow to look at results before "
                   "proceeding (e.g., to look at waveform misfits before "
                   "evaluating the gradient). To determine available options "
                   "for your given workflow: > seisflows print flow"
        )
        self.required.par(
            "SAVEMODEL", required=False, default=True, par_type=bool,
            docstr="Save updated model files after each iteration"
        )
        # Define the Paths required by this module
        self.required.path(
            "FUNC", required=False,
            default=os.path.join(self.path.WORKDIR, "scratch", "evalfunc"),
            docstr="scratch path to store data related to misfit function "
                   "evaluations that take place during the line search. Data "
                   "stored here include residuals from data-synthetic misfit, "
                   "and a given 'try' model being used to generate synthetics."
        )
        # !!! Currently not used
        self.required.path(
            "HESS", required=False,
            default=os.path.join(self.path.WORKDIR, "scratch", "evalhess"),
            docstr="scratch path to store data related to Hessian evaluations"
        )
        self.required.path(
            "OPTIMIZE", required=False,
            default=os.path.join(self.path.WORKDIR, "scratch", "optimize"),
            docstr="scratch path to store data related to nonlinear "
                   "optimization library. Data stored here include model, "
                   "gradient, and search direction vectors (numpy arrays), and"
                   "additional arrays related to specific optimization "
                   "algorithms"
        )

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)

        assert(1 <= self.par.BEGIN <= self.par.END), \
            f"Incorrect BEGIN or END parameter. Values must be in order: " \
            f"1 <= {self.par.BEGIN} <= {self.par.END}"

    def setup(self, flow=None, return_flow=False):
        """
        Lays groundwork for inversion by running setup() functions for the
        involved sub-modules, generating True model synthetic data if necessary,
        and generating the pre-requisite database files.

        .. note::
            This function should only be run one time, at the start of iter 1
        """
        super().setup(flow=flow, return_flow=return_flow)

        self.module("optimize").setup()

    def finalize(self):
        """
        Saves results from current model update iteration and increment the
        iteration number to set up for the next iteration. Finalization is
        expected to the be LAST function in workflow.main()'s  flow list.
        """
        self.logger.info(msg.mjr(f"FINALIZING ITERATION {optimize.iter}"))

        self.checkpoint()
        preprocess.finalize()

    def main(self, flow=None, return_flow=False):
        """
        Overwrites the forward() main function to provide the ability to run
        multiple iterations in a single workflow.

        :type flow: list or tuple
        :param flow: list of Class methods that will be run in the order they
            are provided.
        :type return_flow: bool
        :param return_flow: for CLI tool, simply returns the flow function
            rather than running the workflow. Used for print statements etc.
        """
        if flow is None:
            flow = (self.evaluate_initial_misfit,
                    self.evaluate_gradient,
                    self.process_kernels,
                    self.write_gradient,
                    self.compute_direction,
                    self.line_search,
                    self.export,
                    self.clean
                    )

        self.setup(flow, return_flow)
        optimize = self.module("optimize")

        # Run the workflow until from the current iteration until PAR.END
        optimize.iter = self.par.BEGIN
        self.logger.info(msg.mjr("STARTING INVERSION WORKFLOW"))
        while True:
            self.logger.info(
                msg.mnr(f"ITERATION {optimize.iter} / {self.par.END}")
            )

            # Execute the functions within the flow
            for func in flow[self.start:self.stop]:
                func()
            self.logger.info(msg.mjr(f"FINISHED FLOW EXECUTION"))

            # Reset flow for subsequent iterations
            self.start, self.stop = None, None
            if optimize.iter >= self.par.END:
                break
            optimize.iter += 1

            self.logger.info(msg.sub(f"INCREMENT ITERATION TO {optimize.iter}"))

        self.logger.info(msg.mjr("FINISHED INVERSION WORKFLOW"))

    def evaluate_initial_misfit(self):
        """Inherits from seisflows.workflow.forward.Forward"""
        self.evaluate_initial_misfit()

    def compute_direction(self):
        """
        Computes search direction using the optimization library
        """
        optimize = self.module("optimize")
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
        optimize = self.module("optimize")

        # Calculate the initial step length based on optimization algorithm
        if optimize.line_search.step_count == 0:
            self.logger.info(msg.mjr(f"CONDUCTING LINE SEARCH: "
                                     f"i{optimize.iter:0>2}"
                                     f"s{optimize.line_search.step_count:0>2}")
                             )
            optimize.initialize_search()

        # Attempt a new trial step with the given step length
        optimize.line_search.step_count += 1
        self.logger.info(msg.mnr(f"TRIAL STEP COUNT: "
                                 f"i{optimize.iter:0>2}"
                                 f"s{optimize.line_search.step_count:0>2}"))
        self._evaluate_function(path=self.path.FUNC, suffix="try")

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

    def clean(self):
        """
        Cleans directories in which function and gradient evaluations were
        carried out
        """
        self.logger.info(msg.mnr("CLEANING WORKDIR FOR NEXT ITERATION"))

        unix.rm(self.path.GRAD)
        unix.rm(self.path.FUNC)
        unix.mkdir(self.path.GRAD)
        unix.mkdir(self.path.FUNC)

    def export(self):
        """
        Exports various quantities to PATH.OUTPUT (to disk) from the SCRATCH
        directory as SCRATCH is liable to be overwritten at any point of the
        workflow. This takes place at the end of each iteration, before
        the clean() function is called.
        """
        optimize = self.module("optimize")

        if self.par.SAVEMODEL:
            src = optimize.load("m_new")
            dst = os.path.join(self.path.OUTPUT, f"model_{optimize.iter:04d}")
            self.logger.debug(f"exporting model 'm_new' to disk")
            self._write_vector(src, dst)

        if self.par.SAVEGRADIENT:
            src = optimize.load("g_old")
            dst = os.path.join(self.path.OUTPUT, f"grad_{optimize.iter:04d}")
            self._write_vector(src, dst)

        if self.par.SAVEKERNELS:
            src = os.path.join(self.path.GRAD, "kernels")
            dst = os.path.join(self.path.OUTPUT, f"kernels_{optimize.iter:04d}")
            self.logger.debug(f"saving kernels to path:\n{dst}")
            unix.mv(src, dst)

        if self.par.SAVETRACES:
            self.save_traces()

        if self.par.SAVERESIDUALS:
            src = os.path.join(self.path.GRAD, "residuals")
            dst = os.path.join(self.path.OUTPUT,
                               f"residuals_{optimize.iter:04d}")
            unix.mv(src, dst)

    def _write_vector(self, vector, path):
        """
        Convenience function to write vectors as numpy arrays or as model files
        to a given path
        """
        solver = self.module("solver")

        if self.par.SAVEAS in ["binary", "both"]:
            solver.save(solver.split(vector), path)
        if self.par.SAVEAS in ["vector", "both"]:
            np.save(file=path, arr=vector)
