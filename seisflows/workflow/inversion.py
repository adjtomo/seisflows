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

    def __init__(self, start=1, end=1, thrifty=False, optimize="LBFGS",
                 export_model=True, path_eval_func=None, **kwargs):
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

        # Grab iteration from state file, or set None to have setup() set it
        if "iteration" in self._states:
            self.iteration = int(self._states["iteration"])
        else:
            self.iteration = None

    @property
    def evaluation(self):
        """
        Convenience string return for log messages that gives the iteration
        and step count of the current evaluation as a formatted string
        e.g., i01s00
        """
        return f"i{self.iteration:0>2}s{self.optimize.step_count:0>2}"

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
        return [self.generate_synthetic_data,
                self.evaluate_initial_misfit,
                self.run_adjoint_simulations,
                self.postprocess_event_kernels,
                self.evaluate_gradient_from_kernels,
                self.initialize_line_search,
                self.evaluate_line_search_misfit,
                self.update_line_search,
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

        if self.iteration:
            assert(self.start <= self.iteration <= self.end), (
                f"`workflow.iteration`=={self.iteration} must be between "
                f"parameters `start` and `end`"
            )

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
        unix.mkdir(os.path.join(self.path.eval_grad, "residuals"))
        unix.mkdir(os.path.join(self.path.eval_func, "residuals"))

        self.optimize = self._modules.optimize  # NOQA

        if self.iteration is None:
            self.iteration = self.start

    def run(self):
        """Call the forward.run() function iteratively, from `start` to `end`"""
        while self.iteration < self.end + 1:
            logger.info(msg.mjr(f"RUNNING ITERATION {self.iteration:0>2}"))
            super().run()  # Runs task list
            # Assuming that if `stop_after` is used, that we are NOT iterating
            if self.stop_after is None:
                logger.info(msg.mjr(f"COMPLETE ITERATION {self.iteration:0>2}"))
                self.iteration += 1
                logger.info(f"setting current iteration to: {self.iteration}")
                # Set the state file to pending for new iteration
                self._states = {key: 0 for key in self._states}
                self.checkpoint()
            else:
                break

    def checkpoint(self):
        """
        Add an additional line in the state file to keep track of iteration
        """
        super().checkpoint()

        with open(self.path.state_file, "r") as f:
            lines = f.readlines()
        # Clear out the previous 'iteration' line by reverse scanning through
        # the list and deleting any mention of 'iteration'
        for i in range(len(lines)-1, -1, -1):
            if lines[i].startswith("iteration"):
                lines.pop(i)
        lines.append(f"iteration: {self.iteration}")

        # Rewrite checkpoint file with new iteration line
        with open(self.path.state_file, "w") as f:
            f.writelines(lines)

    def generate_synthetic_data(self, **kwargs):
        """
        Function Override of `workflow.forward.generate_synthetic_data` 

        Add an additional criteria (iteration > 1) that skips over this function
        """
        # We only want to generate synthetic data one time
        if self.iteration > 1:
            logger.debug("inversion iteration > 1, skipping synthetic data gen")
            return
        
        super().generate_synthetic_data(**kwargs)

    def evaluate_objective_function(self, save_residuals=False, components=None,
                                    **kwargs):
        """
        Function Override of `workflow.forward.evaluate_objective_function`

        Simple override to include iteration and step count parameters into
        preprocessing for file naming and tagging. Machinery remains the same.

        .. note::
            Must be run by system.run() so that solvers are assigned individual
            task ids/ working directories.

        :type save_residuals: str
        :param save_residuals: if not None, path to write misfit/residuls to
        :type components: list
        :param components: optional list of components to ignore preprocessing
            traces that do not have matching components. The adjoint sources for
            these components will be 0. E.g., ['Z', 'N']. If None, all available
            components will be considered.
        """
        super().evaluate_objective_function(save_residuals=save_residuals,
                                            components=components,
                                            iteration=self.iteration,
                                            step_count=self.optimize.step_count,
                                            **kwargs
                                            )

    def sum_residuals(self, residuals_files, save_to):
        """
        Convenience function to read in text files containing misfit residual
        information written by `preprocess.quantify_misfit` for each event, and
        sum the total misfit for the evaluation in a given optimization vector.

        Follows Tape et al. 2010 equations 6 and 7

        :type residuals_files: list of str
        :param residuals_files: pathnames to residuals files for each source,
            generated by the preprocessing module. Will be read in and summed
            to provide total misfit
        :type save_to: str
        :param save_to: name of Optimization module vector to save the misfit 
            value 'f', options are 'f_new' for misfit of current accepted model
            'm_new', or 'f_try' for the misfit of the current line search trial
            model
        :rtype: float
        :return: sum of squares of residuals, total misfit
        """
        # Catch empty files because usually we feed this function with a glob
        if not residuals_files:
            logger.critical(
                msg.cli(f"Residuals files not found for {self.evaluation}, "
                        f"preprocessing may have failed. Please check logs.",
                        border="=", header="preprocessing failed")
                )
            sys.exit(-1)

        event_misfits = []
        for residuals_file in residuals_files:
            event_misfit = np.loadtxt(residuals_file)
            # Some preprocessing modules only return a single misfit value
            # which will fail when called with len()
            try:
                num_measurements = len(event_misfit)
            except TypeError:
                num_measurements = 1
            # Tape et al. (2010) Equation 6
            event_misfit = np.sum(event_misfit) / (2. * num_measurements)
            event_misfits.append(event_misfit)

        # Tape et al. (2010) Equation 7
        total_misfit = np.sum(event_misfits) / len(event_misfits)

        # Sum newly calc'd residuals into the optimization library
        self.optimize.save_vector(name=save_to, m=total_misfit)
        logger.info(f"misfit {save_to} ({self.evaluation}) = "
                    f"{total_misfit:.3E}")

    def evaluate_initial_misfit(self, path_model=None, save_residuals=None,
                                **kwargs):
        """
        Overwrite `workflow.forward` to skip over initial misfit evaluation
        (using `MODEL_INIT`) if we are past iteration 1. Additionally, sum
        residuals output by preprocess module and save float to disk, to be
        discoverable by the optimization library

        :type path_model: str
        :param path_model: path to the model files that will be used to evaluate
            initial misfit. If not given, defaults to searching for model
            provided in `path_model_init`.
        :type save_residuals: str
        :param save_residuals: Location to save 'residuals_*.txt files which are
            used to calculate total misfit (f_new).
            - Requires a string formatter '{src}', e.g., 'residual_{src}.txt'
            - String formatter used by preprocessing module to tag files for
            each source to avoid multiple processes writing to the same file.
            - Remainder of string may be some combination of the
            iteration, step count etc. Determined by calling workflow.

        Keyword Arguments
        ::
            bool sum_residuals:
                Bool to determine whether to sum all residuals files saved under
                `save_residuals` filenames. The default behavior should be True,
                that is, once we run preprocessing, we should calculate the
                misfit 'f'. This flag is an option because some workflows, like
                ambient noise inversion, require multiple forward runs before
                summing residuals, so it doesn't make sense to sum each time
                this function is called.
        """
        # Set behavior and paths for saving the residuals_*.txt files
        sum_residuals = kwargs.get("sum_residuals", True)

        # Inversion workflows tag residual file by iteration and step count.
        # Since this is the initial misfit, we assume the step count == 0
        if save_residuals is None:
            save_residuals = os.path.join(
                self.path.eval_grad, "residuals", 
                f"residuals_{{src}}_{self.evaluation}.txt"
            )
        else:
            assert("{src}" in save_residuals), (
                "Workflow path `save_residuals` requires string formatter "
                "{src} within the string name"
            )

        if self.iteration == 1:
            super().evaluate_initial_misfit(save_residuals=save_residuals,
                                            **kwargs)

            # Expose the initial model to the optimization library
            model = Model(self.path.model_init,
                          parameters=self.solver._parameters,
                          regions=self.solver._regions  # req. for 3D_GLOBE only
                          )
            self.optimize.save_vector(name="m_new", m=model)
        else:
            # Thrifty inversion SKIPS initial misfit evaluation, re-using final
            # model from previous line search. This can only happen mid-workflow
            # because we are assuming no solver/preprocess parameters have
            # changed. If they have, then we need to rerun the fwd solver
            if self.thrifty and self._thrifty_status:
                logger.info(msg.mnr("THRIFTY INVERSION; SKIP MISFIT EVAL"))
                return
            # Non-Thrifty, run forward simulation with previous model.
            else:
                logger.info(msg.mnr("EVALUATING MISFIT FOR MODEL `m_new`"))

                # Previous line search will have saved `m_new` as the initial
                # model, export in SPECFEM format to be discoverable by solver
                path_model = os.path.join(self.path.eval_grad, "model")
                m_new = self.optimize.load_vector("m_new")
                m_new.write(path=path_model)

                super().evaluate_initial_misfit(path_model=path_model,
                                                save_residuals=save_residuals,
                                                **kwargs)

        # Read in all avail residual files generated by preprocessing module and
        # sum together to write misfit for initial model: 'f_new'
        if sum_residuals:
            residuals_files = glob(save_residuals.format(src="*"))
            self.sum_residuals(residuals_files, save_to="f_new")
        else:
            logger.info("skipping misfit summation")

    def run_forward_simulations(self, path_model, save_traces=None,
                                export_traces=None, **kwargs):
        """
        Overrides 'workflow.forward.run_forward_simulation' to hijack the
        default path location for exporting traces to disk
        """
        if export_traces is None:
            # e.g., output/i01s00/{source}/syn/*
            export_traces = os.path.join(
                self.path.output, "solver", self.evaluation,
                self.solver.source_name, f"syn"
            )

        super().run_forward_simulations(path_model, save_traces=save_traces,
                                        export_traces=export_traces, **kwargs
                                        )

    def _run_adjoint_simulation_single(self, save_kernels=None,
                                       export_kernels=None, **kwargs):
        """
        Overrides 'workflow.migration._run_adjoint_simulation_single' to hijack
        the default path location for exporting kernels to disk
        """
        # Set default value for `export_kernels` or take program default
        if export_kernels is None:
            export_kernels = os.path.join(
                self.solver.path._solver_output, "kernels", self.evaluation,
                self.solver.source_name
            )

        super()._run_adjoint_simulation_single(save_kernels=save_kernels,
                                               export_kernels=export_kernels,
                                               **kwargs)

    def evaluate_gradient_from_kernels(self):
        """
        Overwrite `workflow.migration` to convert the current model and the
        gradient calculated by migration from their native SPECFEM model format
        into optimization vectors that can be used for model updates.

        Also includes search direction computation, which takes the gradient
        `g_new` and scales to provide an appropriate search direction. At 
        the simplest form (gradient descent), the search direction is simply -g
        """
        super().evaluate_gradient_from_kernels()

        for tag in ["gradient", "kernels"]:
            # e.g., 'gradient' -> 'GRADIENT_01'
            src = os.path.join(self.path.output, tag)
            if os.path.exists(src):
                dst = os.path.join(self.path.output, 
                    f"{tag.upper()}_{self.iteration:0>2}")
                if os.path.exists(dst):
                    unix.rm(dst)
                logger.debug(f"{src} -> {dst}")
                unix.mv(src, dst)

        # Expose the gradient to the optimization library
        logger.info("exposing gradient `g_new` to optimization library")
        gradient = Model(path=os.path.join(self.path.eval_grad, "gradient"),
                         regions=self.solver._regions
                         )
        self.optimize.save_vector(name="g_new", m=gradient)

        # Compute search direction `p_new`: P is used to perturb starting model
        # and is required for the line search 
        logger.info("calculating search direction `p_new` from gradient")
        p_new = self.optimize.compute_direction()

        # Save optimization updates to disk for restarts/checkpoint
        self.optimize.checkpoint()
        self.optimize.save_vector(name="p_new", m=p_new)

    def initialize_line_search(self):
        """
        Computes search direction using the optimization library and sets up
        line search machinery to 'perform line search' by placing correct files
        on disk for each of the modules to find.

        Optimization module perturbs the current model (m_new) by the search
        direction (p_new) to recover the trial model (m_try). This model is
        then exposed on disk to the solver.
        """
        # Set up the line search machinery. Step count forced to 1
        self.optimize.initialize_search()

        # Determine the model we will use for first line search step count
        # m_try = m_new + alpha * p_new (i.e., m_i+1 = m_i + dm)
        alpha, _ = self.optimize.calculate_step_length()
        m_try = self.optimize.compute_trial_model(alpha=alpha)
        
        # Save the current state of the optimization module to disk
        self.optimize.save_vector(name="m_try", m=m_try)
        self.optimize.save_vector(name="alpha", m=alpha)
        self.optimize.checkpoint()

        # Expose model `m_try` to the solvers by placing it in eval_func dir.
        _path_m_try = os.path.join(self.path.eval_func, "model")
        m_try.write(path=_path_m_try)

    def evaluate_line_search_misfit(self):
        """
        Evaluate line search misfit `f_try` by running forward simulations 
        through the trial model `m_try` and comparing with observations.
        Acts like a stripped down version of `evaluate_initial_misfit`

        TODO Add in export traces functionality, need to honor step count
        """
        logger.info(msg.sub(f"LINE SEARCH STEP COUNT "
                            f"{self.optimize.step_count:0>2}"))
        
        logger.info(f"`m_try` model parameters for line search evaluation:")
        self.solver.check_model_values(path=os.path.join(self.path.eval_func, 
                                                         "model"))

        self.system.run(
            [self.run_forward_simulations,
             self.evaluate_objective_function],
            path_model=os.path.join(self.path.eval_func, "model"),
            save_residuals=os.path.join(
                self.path.eval_func, "residuals",
                f"residuals_{{src}}_{self.evaluation}.txt")
        )

        residuals_files = glob(os.path.join(
            self.path.eval_func, "residuals", 
            f"residuals_*_{self.evaluation}.txt")
            )
        assert residuals_files, (
                f"No residuals files found for evaluation {self.evaluation} "
                f"Please check preprocessing log files."
                )
        # Read in all avail residual files generated by preprocessing module
        self.sum_residuals(residuals_files, save_to="f_try")
    
    def update_line_search(self):
        """
        Given the misfit `f_try` calculated in `evaluate_line_search_misfit`,
        use the Optimization module to determine if the line search has passed,
        failed, or needs to perform a subsequent step.

        The line search state machine acts in the following way:
        - Pass: Run clean up and proceed with workflow
        - Try: Re-calculate step length (alpha) and re-evaluate misfit (f_try)
        - Fail: Try to restart optimization module and restart line search. If
            still failing, exit workflow

        .. note::

            Line search starts on step_count == 1 because step_count == 0 is
            considered the misfit of the starting model
        """
        # Update line search history with the step length (alpha) and misfit (f)
        # and incremement the step count
        self.optimize.update_search()
        alpha, status = self.optimize.calculate_step_length()
        m_try = self.optimize.compute_trial_model(alpha=alpha)

        # Save new model (m_try) and step length (alpha) for new trial step
        if alpha is not None:
            self.optimize.save_vector("alpha", alpha)
        if m_try is not None:
            self.optimize.save_vector("m_try", m_try)

        # Proceed based on the outcome of the line search
        if status.upper() == "PASS":
            # Save outcome of line search to disk; reset step to 0 for next iter
            logger.info("trial step successful. finalizing line search")

            # Finalizing line search sets `m_try` -> `m_new` for later iters
            self.optimize.finalize_search()
            self.optimize.checkpoint()
            return
        elif status.upper() == "TRY":
            logger.info("trial step unsuccessful. re-attempting line search")

            # Expose the new model to the solver directories for the next step
            _path_m_try = os.path.join(self.path.eval_func, "model")
            m_try.write(path=_path_m_try)

            # Re-set state file to ensure that job failure will recover
            self._states["evaluate_line_search_misfit"] = 0

            # Recursively run the line search to get a new misfit
            self.optimize.checkpoint()
            self.evaluate_line_search_misfit()
            self.update_line_search()  # RECURSIVE CALL
        elif status.upper() == "FAIL":
            # Check if we are able to restart line search w/ new parameters
            if self.optimize.attempt_line_search_restart():
                logger.info("line search has failed. restarting "
                            "optimization algorithm and line search.")
                # Reset the line search machinery; set step count to 0
                self.optimize.restart()

                # Re-set state file to ensure that job failure will recover
                self._states["evaluate_line_search_misfit"] = 0

                # Restart the entire line search procedure
                self.optimize.checkpoint()
                self.initialize_line_search()
                self.evaluate_line_search_misfit()
                self.update_line_search()  # RECURSIVE CALL
            # If we can't then line search has failed. Abort workflow
            else:
                logger.critical(
                    msg.cli("Line search has failed to reduce the misfit and "
                            "has run out of fallback options. Aborting "
                            "inversion.", border="=",
                            header="line search failed")
                )
                sys.exit(-1)

    def finalize_iteration(self):
        """
        Cleans directories in which function and gradient evaluations were
        carried out. Contains some logic to consider whether or not to continue
        with a thrifty inversion.
        """
        super().finalize_iteration()

        # Export scratch files to output if requested
        if self.export_model:
            model = self.optimize.load_vector("m_new")
            m_fid = os.path.join(self.path.output, 
                                 f"MODEL_{self.iteration:0>2}")
            logger.info(f"writing model `m_new` to {m_fid}")
            model.write(path=m_fid)

        # Organize log files to keep file count in the main log dir. low
        logger.debug(f"organizing log files in {self.system.path.log_files}")
        logdir = os.path.join(self.system.path.log_files, 
                              f"LOGS_i{self.iteration:0>2}")
        unix.mkdir(logdir)
        fids = glob(os.path.join(self.system.path.log_files, "*"))
        fids = [fid for fid in fids if not 
                os.path.basename(fid).startswith("LOGS_")]
        unix.mv(fids, logdir)

        # Update optimization on disk 
        self.optimize.checkpoint()

        # Thrifty Inversion keeps last function evaluation for next iteration
        self._thrifty_status = self._update_thrifty_status()
        if self._thrifty_status:
            unix.rm(self.path.eval_grad)
            # Eval func model now defines the current model 'm_new'
            unix.mv(self.path.eval_func, self.path.eval_grad)
            unix.mkdir(self.path.eval_func)
            unix.mkdir(os.path.join(self.path.eval_func, "residuals"))
        # Std. Inversion bombs out both eval dirs., next iter will start fresh
        else:
            logger.info("cleaning out scratch directory for next iteration")
            unix.rm(self.path.eval_grad)
            unix.rm(self.path.eval_func)

            unix.mkdir(self.path.eval_grad)
            unix.mkdir(self.path.eval_func)

            unix.mkdir(os.path.join(self.path.eval_grad, "residuals"))
            unix.mkdir(os.path.join(self.path.eval_func, "residuals"))

    def _update_thrifty_status(self):
        """
        Determine if line search forward simulation can be carried over to the
        next iteration. Checks criteria related to the current iteration and
        its position relative to the start and end of the workflow.

        .. note::
            Resumed, failed workflows will not re-load `_thrifty_status` so
            initial misfit will always be evaluated in that case.
        """
        # !!! Forcing first iteration to be thrifty because i'm not going to
        # !!! change the parameters
        if self.iteration == self.start:
            _thrifty_status = True
        # if self.iteration == self.start:
        #     logger.info("thrifty inversion encountering first iteration, "
        #                 "defaulting to standard inversion workflow")
        #     _thrifty_status = False
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
