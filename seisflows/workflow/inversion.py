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

from seisflows import logger
from seisflows.workflow.migration import Migration
from seisflows.tools import msg, unix
from seisflows.tools.specfem import Model


class Inversion(Migration):
    """
    [workflow.inversion] Peforms iterative nonlinear inversion using a built-in
    optimization library which stores model vectors on disk.

    :type start: int
    :param start: start inversion workflow at this iteration. 1 <= start <= inf
    :type end: int
    :param end: end inversion workflow at this iteration. start <= end <= inf
    :type export_model: bool
    :param export_model: export best-fitting model from the line search to disk.
        If False, new models can be discarded from scratch at any time.
    :type path_eval_func: str
    :param path_eval_func: scratch path to store files for line search objective
        function evaluations, including models, misfit and residuals
    """
    __doc__ = Migration.__doc__ + __doc__

    def __init__(self, modules=None, start=1, end=1, export_model=True,
                 path_eval_func=None, **kwargs):
        """Instantiate Inversion-specific parameters"""

        super().__init__(**kwargs)

        self._modules = modules
        self.start = start
        self.end = end
        self.export_model = export_model

        self.path.eval_func = path_eval_func or \
                              os.path.join(self.path.workdir, "scratch",
                                           "eval_func")

        # Overwriting base class required modules list
        self._required_modules = ["system", "solver", "preprocess", "optimize"]

        # Empty module variables that should be filled in by setup
        self.optimize = None

    @property
    def task_list(self):
        """
        USER-DEFINED TASK LIST. This property defines a list of class methods
        that take NO INPUT and have NO RETURN STATEMENTS. This defines your
        linear workflow, i.e., these tasks are to be run in order from start to
        finish to complete a workflow.

        This excludes 'check' (which is run during 'import_seisflows') and
        'setup' which should be run separately

        .. note::
            For workflows that require an iterative approach (e.g. inversion),
            this task list will be looped over, so ensure that any setup and
            teardown tasks (run once per workflow, not once per iteration) are
            not included.

        :rtype: list
        :return: list of methods to call in order during a workflow
        """
        return [self.evaluate_initial_misfit,
                self.run_adjoint_simulations,
                self.postprocess_event_kernels,
                self.evaluate_gradient_from_kernels,
                self.compute_search_direction,
                self.perform_line_search,
                self.clean_scratch_directory
                ]

    def check(self):
        """
        Checks inversion-specific parameters
        """
        super().check()

        assert(1 <= self.start <= self.end), \
            f"Incorrect START or END parameter. Values must be in order: " \
            f"1 <= {self.start} <= {self.end}"

    def setup(self):
        """
        Lays groundwork for inversion by running setup() functions for the
        involved sub-modules, generating True model synthetic data if necessary,
        and generating the pre-requisite database files.
        """
        super().setup()

        unix.mkdir(self.path.eval_func)
        self.optimize = self._modules.optimize

    def run(self):
        """Call the forward.run() function iteratively, from `start` to `end`"""
        for i in range(self.start, self.end + 1):
            logger.info(msg.mnr(f"RUNNING ITERATION {i:0>2}"))
            super().run()
            logger.info(msg.mnr(f"COMPLETED ITERATION {i:0>2}"))

    def evaluate_initial_misfit(self):
        """
        Overwrite `workflow.forward` just to sum residuals output by preprocess
        module and save them to disk, to be discoverable by the optimization
        library
        """
        super().evaluate_initial_misfit()

        # Override function to sum residuals into the optimization library
        residuals = np.loadtxt(os.path.join(self.path.eval_grad, "residuals"))
        total_misfit = self.preprocess.sum_residuals(residuals)
        self.optimize.save(name="f_new", m=total_misfit)

    def evaluate_gradient_from_kernels(self):
        """
        Overwrite `workflow.migration` to convert the current model and the
        gradient calculated by migration from their native SPECFEM model format
        into optimization vectors that can be used for model updates.
        """
        super().evaluate_gradient_from_kernels()

        model = Model(os.path.join(self.path.eval_grad, "model"))
        self.optimize.save(name="m_new", m=model)

        gradient = Model(path=os.path.join(self.path.eval_grad, "gradient"))
        self.optimize.save(name="g_new", m=gradient)

    def compute_search_direction(self):
        """
        Computes search direction using the optimization library and performs
        a sets up line search machinery to 'perform line search'
        """
        logger.info(f"computing search direction with "
                    f"'{self.optimize.__class__.__name__}'")
        p_new = self.optimize.compute_direction()
        self.optimize.save(name="p_new", m=p_new)

        # Check that our search direction will actually perturb the model
        if sum(p_new.vector) == 0:
            logger.critical(msg.cli(
                "Search direction vector 'p' is 0, meaning no model update can "
                "take place. Please check your gradient and waveform misfits. "
                "SeisFlows exiting prior to start of line search.", border="=",
                header="line search error")
            )
            sys.exit(-1)

        logger.info(
            msg.mnr(f"INITALIZING LINE SEARCH: i{self.optimize.iteration:0>2}"
                    f"s{self.optimize.step_count:0>2}")
                    )
        m_try, alpha = self.optimize.initialize_search()
        self.optimize.save(name="m_try", m=m_try)
        self.optimize.save(name="alpha", m=alpha)

        # Expose model `m_try` to the solver by placing it in eval_func dir.
        m_try.write(path=os.path.join(self.path.eval_func, "model"))

    def perform_line_search(self):
        """
        Conducts line search in given search direction until the objective
        function is reduced acceptably, or the line search fails due to
        user-defined limit criteria.

        Status codes:
            status > 0  : finished
            status == 0 : not finished
            status < 0  : failed
        """
        while True:
            # Attempt a new trial step with the given step length
            self.optimize.line_search.step_count += 1
            logger.info(
                msg.sub(f"TRIAL STEP COUNT: i{self.optimize.iteration:0>2}"
                        f"s{self.optimize.step_count:0>2}")
            )

            # Run the forward simulation on system and calculate misfit as f_try
            self.system.run(
                [self.evaluate_objective_function],
                path_model=os.path.join(self.path.eval_func, "model"),
                save_residuals=os.path.join(self.path.eval_func, "residuals")
            )
            residuals = np.loadtxt(os.path.join(self.path.eval_func,
                                                "residuals"))
            total_misfit = self.preprocess.sum_residuals(residuals)
            self.optimize.save(name="f_try", m=total_misfit)

            # Check the function evaluation against line search history
            status = self.optimize.update_search()

            # Proceed based on the outcome of the line search
            if status > 0:
                # Save outcome of line search to disk; reset step to 0 for next iter
                logger.info("trial step successful. finalizing line search")
                self.optimize.finalize_search()
                return
            elif status == 0:
                logger.info("trial step unsuccessful. retrying with new trial "
                            "step")
                # Recursively call this function to attempt another trial step
                continue
            elif status < 0:
                if self.optimize.retry_status():
                    logger.info("line search has failed. restarting "
                                "optimization algorithm and line search.")
                    # Reset the line search machinery; set step count to 0
                    self.optimize.restart()
                    continue
                else:
                    logger.info("line search has failed. aborting inversion.")
                    sys.exit(-1)

    def clean_scratch_directory(self):
        """
        Cleans directories in which function and gradient evaluations were
        carried out
        """
        logger.info(msg.mnr("CLEANING WORKDIR FOR NEXT ITERATION"))

        unix.rm(self.path.eval_grad)
        unix.rm(self.path.eval_func)

        unix.mkdir(self.path.eval_grad)
        unix.mkdir(self.path.eval_func)

    # def export(self):
    #     """
    #     Exports various quantities to PATH.OUTPUT (to disk) from the SCRATCH
    #     directory as SCRATCH is liable to be overwritten at any point of the
    #     workflow. This takes place at the end of each iteration, before
    #     the clean() function is called.
    #     """
    #     optimize = self.module("optimize")
    #
    #     if self.par.SAVEMODEL:
    #         src = os.path.join(self.path.OPTIMIZE, "m_new")
    #         dst = os.path.join(self.path.OUTPUT, f"model_{optimize.iter:04d}")
    #         logger.debug(f"exporting model 'm_new' to disk")
    #         unix.cp(src, dst)
    #
    #     if self.par.SAVEGRADIENT:
    #         src = os.path.join(self.path.OPTIMIZE, "g_old")
    #         dst = os.path.join(self.path.OUTPUT, f"grad_{optimize.iter:04d}")
    #         unix.cp(src, dst)
    #
    #     if self.par.SAVEKERNELS:
    #         src = os.path.join(self.path.GRAD, "kernels")
    #         dst = os.path.join(self.path.OUTPUT, f"kernels_{optimize.iter:04d}")
    #         logger.debug(f"saving kernels to path:\n{dst}")
    #         unix.mv(src, dst)
    #
    #     if self.par.SAVERESIDUALS:
    #         src = os.path.join(self.path.GRAD, "residuals")
    #         dst = os.path.join(self.path.OUTPUT,
    #                            f"residuals_{optimize.iter:04d}")
    #         unix.mv(src, dst)

