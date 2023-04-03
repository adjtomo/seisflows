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

from glob import glob
from seisflows import logger
from seisflows.workflow.migration import Migration
from seisflows.tools import msg, unix
from seisflows.tools.model import Model


class Inversion(Migration):
    """
    Inversion Workflow
    ------------------
    Peforms iterative nonlinear inversion using the machinery of the Forward
    and Migration workflows, as well as a built-in optimization library.

    Parameters
    ----------
    :type start: int
    :param start: start inversion workflow at this iteration. 1 <= start <= inf
    :type end: int
    :param end: end inversion workflow at this iteration. start <= end <= inf
    :type iteration: int
    :param iteration: The current iteration of the workflow. If NoneType, takes
        the value of `start` (i.e., first iteration of the workflow). User can
        also set between `start` and `end` to resume a failed workflow.
    :type thrifty: bool
    :param thrifty: a thrifty inversion skips the costly intialization step
        (i.e., forward simulations and misfit quantification) if the final
        forward simulations from the previous iterations line search can be
        used in the current one. Requires L-BFGS optimization.
    :type export_model: bool
    :param export_model: export best-fitting model from the line search to disk.
        If False, new models can be discarded from scratch at any time.

    Paths
    -----
    :type path_eval_func: str
    :param path_eval_func: scratch path to store files for line search objective
        function evaluations, including models, misfit and residuals
    ***
    """
    __doc__ = Migration.__doc__ + __doc__

    def __init__(self, modules=None, start=1, end=1,
                 thrifty=False, optimize="LBFGS", export_model=True,
                 path_eval_func=None,
                 **kwargs):
        """
        Instantiate Inversion-specific parameters. Non-essential parameters are
        listed here, rather than in the class docstring.

        :type optimize: str
        :param optimize: Name of the optimization module chosen by the user.
            This should be instantiated by default when using `import_seisflows`
            Used to check that the correct module is set when performing a
            `thrifty` inversion.
        """
        super().__init__(**kwargs)

        self._modules = modules
        self.start = start
        self.end = end
        self.export_model = export_model
        self.thrifty = thrifty

        # Append an additional path for line search function evaluations
        self.path["eval_func"] = path_eval_func or \
                                 os.path.join(self.path.workdir, "scratch",
                                           "eval_func")

        # Internal attribute for keeping track of inversion
        self._optimize_name = optimize
        self._thrifty_status = False
        self._required_modules = ["system", "solver", "preprocess", "optimize"]

        # Grab iteration from state file
        if "iteration" in self._states:
            self.iteration = int(self._states["iteration"])
            logger.debug(f"setting iteration=={self.iteration} from state file")
        else:
            self.iteration = start

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
                self.finalize_iteration
                ]

    def check(self):
        """
        Checks inversion-specific parameters
        """
        super().check()

        assert(1 <= self.start <= self.end), \
            f"Incorrect START or END parameter. Values must be in order: " \
            f"1 <= {self.start} <= {self.end}"

        assert(self.start <= self.iteration <= self.end), \
            f"`workflow.iteration` must be between `start` and `end`"

        if self.iteration > 1:
            assert(os.path.exists(self.path.eval_grad)), \
                f"scratch path `eval_grad` does not exist but should for a " \
                f"workflow with `iteration` >= 1"

        if self.iteration >= self.end + 1:
            logger.warning(f"current `iteration` is >= chosen `end` point. "
                           f"Inversion workflow will not `run`")

        if self.thrifty:
            assert(self._optimize_name == "LBFGS"), (
                f"a `thrifty` inversion requires the optimization module to be "
                f"set as 'LBFGS'"
            )

    def setup(self):
        """
        Assigns modules as attributes of the workflow. I.e., `self.solver` to
        access the solver module (or `workflow.solver` from outside class)

        Lays groundwork for inversion by running setup() functions for the
        involved sub-modules, generating True model synthetic data if necessary,
        and generating the pre-requisite database files.
        """
        super().setup()

        unix.mkdir(self.path.eval_func)

        self.optimize = self._modules.optimize  # NOQA
        # If optimization has been run before, re-load from checkpoint
        self.optimize.load_checkpoint()

    def run(self):
        """Call the forward.run() function iteratively, from `start` to `end`"""
        while self.iteration < self.end + 1:
            logger.info(msg.mnr(f"RUNNING ITERATION {self.iteration:0>2}"))
            super().run()  # Runs task list
            # Assuming that if `stop_after` is used, that we are NOT iterating
            if self.stop_after is None:
                logger.info(msg.mnr(f"COMPLETE ITERATION {self.iteration:0>2}"))
                self.iteration += 1
                logger.info(f"setting current iteration to: {self.iteration}")
                # Set the state file to pending for new iteration
                self._states = {key: "pending" for key in self._states}
                self.checkpoint()
            else:
                break

    def checkpoint(self):
        """
        Add an additional line in the state file to keep track of iteration,
        """
        super().checkpoint()
        with open(self.path.state_file, "r") as f:
            lines = f.readlines()
        # Clear out the previous 'iteration' line and add in new
        for i, line in enumerate(lines[:]):
            if "iteration:" in line:
                lines.pop(i)
                break
        lines.append(f"iteration: {self.iteration}")

        # Rewrite checkpoint file with new iteration line
        with open(self.path.state_file, "w") as f:
            f.writelines(lines)

    def evaluate_objective_function(self, save_residuals=False, **kwargs):
        """
        Overwrite evaluate objective function to include MORE input parameters
        specifying which evaluation in the inversion we are at. Also removes
        the check for a preprocessing module because it is assumed we have a
        preprocsesing module for an inversion workflow.

        .. note::
            Must be run by system.run() so that solvers are assigned individual
            task ids/ working directories.
        """
        logger.debug(f"quantifying misfit with "
                     f"'{self.preprocess.__class__.__name__}'")

        # If line search, add step count as suffix in the residuals file
        if self.path.eval_func == os.path.dirname(save_residuals):
            save_residuals = save_residuals.format(src=self.solver.source_name,
                                                   it=str(self.iteration),
                                                   sc=str(self.optimize.step_count + 1))
        else:
            save_residuals = save_residuals.format(src=self.solver.source_name,
                                                   it=str(self.iteration))

        if self.export_residuals:
            export_residuals = os.path.join(self.path.output, "residuals")
        else:
            export_residuals = False

        self.preprocess.quantify_misfit(
            source_name=self.solver.source_name,
            save_adjsrcs=os.path.join(self.solver.cwd, "traces", "adj"),
            save_residuals=save_residuals,
            export_residuals=export_residuals,
            iteration=self.iteration,
            step_count=self.optimize.step_count,
        )

    def evaluate_initial_misfit(self):
        """
        Overwrite `workflow.forward` to skip over initial misfit evaluation
        (using `MODEL_INIT`) if we are past iteration 1. Additionally, sum
        residuals output by preprocess module and save float to disk, to be
        discoverable by the optimization library
        """
        if self.iteration == 1:
            super().evaluate_initial_misfit()

            # Expose the initial model to the optimization library
            model = Model(self.path.model_init,
                          parameters=self.solver._parameters,
                          regions=self.solver._regions  # 3DGLOBE only
                          )
            self.optimize.save_vector(name="m_new", m=model)
        else:
            # Thrifty inversion SKIPS initial misfit evaluation, re-using final
            # model from previous line search. Can only happen mid-workflow
            if self.thrifty and (
                    self._thrifty_status or self.iteration == self.start):
                logger.info(msg.mnr("THRIFTY INVERSION; SKIP MISFIT EVAL"))
            else:
                logger.info(msg.mnr("EVALUATING MISFIT FOR MODEL `m_new`"))
                # Previous line search will have saved `m_new` as the initial
                # model, export in SPECFEM format to a path discoverable by all
                # solvers
                path_model = os.path.join(self.path.eval_grad, "model")
                m_new = self.optimize.load_vector("m_new")
                m_new.write(path=path_model)

                # Run forward simulation/misfit quantification with previous
                # model
                self.system.run(
                    [self.run_forward_simulations,
                     self.evaluate_objective_function],
                    path_model=path_model,
                    save_residuals=os.path.join(self.path.eval_grad,
                                                "residuals_{src}_{it}.txt")
                )

        # Rename exported synthetic traces so they are not overwritten by
        # future forward simulations
        if self.export_traces:
            unix.mv(src=os.path.join(self.path.output, "solver"),
                    dst=os.path.join(self.path.output,
                                     f"solver_{self.iteration:0>2}"))

        # Override function to sum residuals into the optimization library
        residuals_files = glob(os.path.join(self.path.eval_grad,
                            f"residuals_*_{self.iteration}.txt"))

        residuals = self.preprocess.read_residuals(residuals_files)
        total_misfit = self.preprocess.sum_residuals(residuals)
        self.optimize.save_vector(name="f_new", m=total_misfit)

    def evaluate_gradient_from_kernels(self):
        """
        Overwrite `workflow.migration` to convert the current model and the
        gradient calculated by migration from their native SPECFEM model format
        into optimization vectors that can be used for model updates.
        """
        super().evaluate_gradient_from_kernels()

        # Rename kernels (K) and gradient (G) output files by iteration number
        # so they don't get overwritten by future iterations.
        src = os.path.join(self.path.output, "kernels")
        dst = os.path.join(self.path.output, f"KERNELS_{self.iteration:0>2}")
        if os.path.exists(src):
            logger.debug(f"{src} -> {dst}")
            unix.mv(src, dst)

        src = os.path.join(self.path.output, "gradient")
        dst = os.path.join(self.path.output, f"GRADIENT_{self.iteration:0>2}")
        if os.path.exists(src):
            if os.path.exists(dst):
                unix.rm(dst)
            logger.debug(f"{src} -> {dst}")
            unix.mv(src, dst)

        # Expose the gradient to the optimization library
        gradient = Model(path=os.path.join(self.path.eval_grad, "gradient"),
                         regions=self.solver._regions
                         )
        self.optimize.save_vector(name="g_new", m=gradient)

    def initialize_line_search(self):
        """
        Computes search direction using the optimization library and sets up
        line search machinery to 'perform line search' by placing correct files
        on disk for each of the modules to find.

        Optimization module perturbs the current model (m_new) by the search
        direction (p_new) to recover the trial model (m_try). This model is
        then exposed on disk to the solver.
        """
        logger.info(msg.mnr("RUNNING LINE SEARCH"))
        logger.info(f"initializing "
                    f"'{self.optimize.line_search_method}'ing "
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
        self.optimize.save_vector(name="p_new", m=p_new)

        # Scale search direction with step length alpha generate a model update
        m_try, alpha = self.optimize.initialize_search()
        self.optimize.save_vector(name="m_try", m=m_try)
        self.optimize.save_vector(name="alpha", m=alpha)
        self.optimize.checkpoint()

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
        logger.info(msg.sub(f"LINE SEARCH STEP COUNT "
                            f"{self.optimize.step_count + 1:0>2}"))

        # Run fwd solver with the model 'm_try'. Corresponding misfit is 'f_try'
        self._evaluate_line_search_misfit()

        # Increment step count, calculate new step length/model, check misfit
        m_try, alpha, status = self.optimize.update_line_search()
        self.optimize.checkpoint()

        # Proceed based on the outcome of the line search
        if status.upper() == "PASS":
            # Save outcome of line search to disk; reset step to 0 for next iter
            logger.info("trial step successful. finalizing line search")

            # Save new model (m_try) and step length (alpha) for records
            self.optimize.save_vector("alpha", alpha)
            self.optimize.save_vector("m_try", m_try)
            m_try.write(path=os.path.join(self.path.eval_func, "model"))
            del m_try  # clear potentially large model vector from memory

            self.optimize.finalize_search()
            self.optimize.checkpoint()
            return
        elif status.upper() == "TRY":
            logger.info("trial step unsuccessful. re-attempting line search")

            # Save new model (m_try) and step length (alpha) for new trial step
            self.optimize.save_vector("alpha", alpha)
            self.optimize.save_vector("m_try", m_try)
            m_try.write(path=os.path.join(self.path.eval_func, "model"))
            del m_try  # clear potentially large model vector from memory

            # Checkpoint and re-run line search evaluation
            self.optimize.checkpoint()
            self.perform_line_search()  # RECURSIVE CALL
        elif status.upper() == "FAIL":
            # Check if we are able to restart line search w/ new parameters
            if self.optimize.attempt_line_search_restart():
                logger.info("line search has failed. restarting "
                            "optimization algorithm and line search.")
                # Reset the line search machinery; set step count to 0
                self.optimize.restart()
                self.optimize.checkpoint()
                self.perform_line_search()  # RECURSIVE CALL
            # If we can't then line search has failed. Abort workflow
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
            [self.run_forward_simulations,
             self.evaluate_objective_function],
            path_model=os.path.join(self.path.eval_func, "model"),
            save_residuals=os.path.join(self.path.eval_func,
                                        "residuals_{src}_{it}_{sc}.txt")
        )

        residuals_files = glob(os.path.join(self.path.eval_func,
                            f"residuals_*_{self.iteration}_{self.optimize.step_count + 1}.txt"))

        residuals = self.preprocess.read_residuals(residuals_files)
        total_misfit = self.preprocess.sum_residuals(residuals)
        logger.debug(f"misfit for trial model (f_try) == {total_misfit:.2E}")
        self.optimize.save_vector(name="f_try", m=total_misfit)

    def finalize_iteration(self):
        """
        Cleans directories in which function and gradient evaluations were
        carried out. Contains some logic to consider whether or not to continue
        with a thrifty inversion.
        """
        logger.info(msg.sub("CLEANING WORKDIR FOR NEXT ITERATION"))

        # Export scratch files to output if requested
        if self.export_model:
            model = self.optimize.load_vector("m_new")
            model.write(path=os.path.join(self.path.output,
                                          f"MODEL_{self.iteration:0>2}"),
                        )

        # Update optimization
        self.optimize.checkpoint()

        # Clear out the scratch directory
        self._thrifty_status = self._update_thrifty_status()
        if self._thrifty_status:
            unix.rm(self.path.eval_grad)
            # Eval func model now defines the current model 'm_new'
            unix.mv(self.path.eval_func, self.path.eval_grad)
            unix.mkdir(self.path.eval_func)
        else:
            unix.rm(self.path.eval_grad)
            unix.rm(self.path.eval_func)

            unix.mkdir(self.path.eval_grad)
            unix.mkdir(self.path.eval_func)

        self.preprocess.finalize()

    def _update_thrifty_status(self):
        """
        Determine if line search forward simulation can be carried over to the
        next iteration. Checks criteria related to the current iteration and
        its position relative to the start and end of the workflow.

        .. note::
            Resumed, failed workflows will not re-load `_thrifty_status` so
            initial misfit will always be evaluated in that case.
        """
        if self.iteration == self.start:
            logger.info("thrifty inversion encountering first iteration, "
                        "defaulting to standard inversion workflow")
            _thrifty_status = False
        elif self.optimize._restarted:  # NOQA
            logger.info("optimization has been restarted, defaulting to "
                        "standard inversion workflow")
            _thrifty_status = False
        elif self.iteration == self.end:
            logger.info("thrifty inversion encountering final iteration, "
                        "defaulting to inversion workflow")
            _thrifty_status = False
        else:
            logger.info("acceptable conditions for thrifty inverison, "
                        "continuing with thrifty inversion")
            _thrifty_status = True

        return _thrifty_status
