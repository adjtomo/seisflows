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
                self.initialize_line_search,
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

    def checkpoint(self):
        """
        Override checkpoint and add the optimization iteration parameter
        """
        super().checkpoint()
        with open(self.path.state_file, "a") as f:
            f.write(f"iteration: {self.optimize.iteration}")

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

    def initialize_line_search(self):
        """
        Computes search direction using the optimization library and sets up
        line search machinery to 'perform line search' by placing correct files
        on disk for each of the modules to find.

        Optimization module perturbs the current model (m_new) by the search
        direction (p_new) to recover the trial model (m_try). This model is
        then exposed on disk to the solver.
        """
        logger.info(f"initializing "
                    f"'{self.optimize.line_search.__class__.__name__}'ing "
                    f"line search")

        # 'p' is the search direction used to perturb the initial model
        p_new = self.optimize.compute_direction()
        if sum(p_new.vector) == 0:
            logger.critical(msg.cli(
                "Search direction vector 'p' is 0, meaning no model update can "
                "take place. Please check your gradient and waveform misfits. "
                "SeisFlows exiting prior to start of line search.", border="=",
                header="line search error")
            )
            sys.exit(-1)
        self.optimize.save(name="p_new", m=p_new)

        # Scale search direction with step length alpha generate a model update
        m_try, alpha = self.optimize.initialize_search()
        self.optimize.save(name="m_try", m=m_try)
        self.optimize.save(name="alpha", m=alpha)
        self.optimize.line_search.save_search_history()

        # Expose model `m_try` to the solver by placing it in eval_func dir.
        m_try.write(path=os.path.join(self.path.eval_func, "model"))

    def perform_line_search(self):
        """
        Conducts line search in given search direction until the objective
        function is reduced acceptably, or the line search fails due to
        user-defined limit criteria.

        .. note::
            Starts on step_count == 1 because step_count == 0 will be the
            misfit of the starting model

        Status codes:
            status > 0  : finished
            status == 0 : not finished
            status < 0  : failed
        """
        # self.optimize.line_search.load_search_history()

        logger.info(msg.sub(f"LINE SEARCH STEP COUNT "
                            f"{self.optimize.step_count + 1:0>2}"))

        # Run fwd solver with the model 'm_try'. Corresponding misfit is 'f_try'
        self._evaluate_line_search_misfit()

        # Increment step count, calculate new step length/model, check misfit
        status = self.optimize.update_line_search()
        self.optimize.line_search.save_search_history()

        # Proceed based on the outcome of the line search
        if status == 1:
            # Save outcome of line search to disk; reset step to 0 for next iter
            logger.info("trial step successful. finalizing line search")
            self.optimize.finalize_search()
            self.optimize.line_search.save_search_history()
            return
        elif status == 0:
            logger.info("trial step unsuccessful. re-attempting line search")
            self.perform_line_search()  # RECURSIVE CALL
        elif status == -1:
            if self.optimize.attempt_line_search_restart():
                logger.info("line search has failed. restarting "
                            "optimization algorithm and line search.")
                # Reset the line search machinery; set step count to 0
                self.optimize.restart()
                self.perform_line_search()  # RECURSIVE CALL
            else:
                logger.critical(
                    msg.cli("Line search has failed to reduce the misfit and "
                            "has run out of fallback options. Aborting "
                            "inversion.", border="=",
                            header="line search failed")
                )
                sys.exit(-1)

    def _evaluate_line_search_misfit(self):
        """Convenience fuinction to wrap forward solver and misfit calc"""
        self.system.run(
            [self.evaluate_objective_function],
            path_model=os.path.join(self.path.eval_func, "model"),
            save_residuals=os.path.join(self.path.eval_func, "residuals")
        )
        residuals = np.loadtxt(os.path.join(self.path.eval_func,
                                            "residuals"))
        total_misfit = self.preprocess.sum_residuals(residuals)
        logger.debug(f"misfit for trial model (f_try) == {total_misfit:.2E}")
        self.optimize.save(name="f_try", m=total_misfit)

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


