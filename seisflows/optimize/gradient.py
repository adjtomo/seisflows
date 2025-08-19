#!/usr/bin/env python3
"""
Gradient descent nonlinear optimization algorithm. Acts as the Base class for
optimization.

The Optimization library contains classes and methods used to solve nonlinear
optimization problems, i.e., misfit minimization. Various subclasses implement
different optimization algorithms.

.. note::
    To reduce memory overhead, and to enable optimization restarts due to
    failed line searches, model/gradient vectors and line search history are
    passed on disk rather than in memory.

.. note::
    The default numerical parameters for each algorithm should work well for a
    range of applications without manual tuning. If the nonlinear
    optimization procedure stagnates, it may be due to issues involving:

    1) poor data quality
    2) choice of objective function
    3) data processing parameters (i.e., filter bandpass)
    4) regularization methods

    Problems in any of these areas usually manifest themselves through
    stagnation of the nonlinear optimizationalgorithm.
"""
import os
import sys
import numpy as np
from glob import glob

from seisflows import logger
from seisflows.tools import msg, unix
from seisflows.tools.config import Dict
from seisflows.tools.graphics import plot_optim_stats
from seisflows.tools.specfem_model import Model
from seisflows.plugins import line_search as line_search_dir


class Gradient:
    """
    Gradient Optimization [Optimize Base]
    -------------------------------------
    Defines foundational structure for Optimization module. Applies a 
    gradient/steepest descent optimization algorithm.

    Parameters
    ----------
    :type line_search_method: str
    :param line_search_method: chosen line_search algorithm. Currently available
        are 'bracket' and 'backtrack'. See seisflows.plugins.line_search
        for all available options
    :type preconditioner: str
    :param preconditioner: algorithm for preconditioning gradients. Currently
        available: 'diagonal'. Requires `path_preconditioner` to point to a
        set of files that define the preconditioner, formatted the same as the
        input model
    :type step_count_max: int
    :param step_count_max: maximum number of trial steps to perform during
        the line search before a change in line search behavior is
        considered, or a line search is considered to have failed.
    :type step_len_init: float
    :param step_len_init: initial line search step length guess, provided
        as a fraction of current model parameters.
    :type step_len_max: float
    :param step_len_max: optional, maximum allowable step length during the line
        search. Set as a fraction of the current model parameters
    :type step_len_min: float
    :param step_len_min: optional, minimum allowable step length during the line
        search. Set as a fraction of the current model parameters

    Paths
    -----
    :type path_preconditioner: str
    :param path_preconditioner: optional path to a set of preconditioner files
        formatted the same as the input model (or output model of solver).
        Required to exist and contain files if `preconditioner`==True
    ***
    """
    def __init__(self, line_search_method="bracket",
                 preconditioner=None, step_count_max=10, step_len_init=0.01,
                 step_len_max=0.1, step_len_min=1E-3, workdir=os.getcwd(), 
                 path_optimize=None, path_output=None, path_preconditioner=None, 
                 **kwargs):
        """
        Gradient-descent input parameters.

        .. note::
            Paths listed here are shared with `workflow.forward` and so are not
            included in the class docstring.

        :type workdir: str
        :param workdir: working directory in which to look for data and store
            results. Defaults to current working directory
        :type path_output: str
        :param path_output: path to directory used for permanent storage on disk.
            Results and exported scratch files are saved here.
        """
        super().__init__()

        self.preconditioner = preconditioner
        self.step_count_max = step_count_max
        self.step_len_init = step_len_init
        self.step_len_max = step_len_max
        self.step_len_min = step_len_min

        # Set required path structure
        self.path = Dict(
            scratch=path_optimize or
                    os.path.join(workdir, "scratch", "optimize"),
            output=path_output or os.path.join(workdir, "output"),
            preconditioner=path_preconditioner,
        )

        # Hidden paths to store checkpoint file in scratch directory
        self.path["_checkpoint"] = os.path.join(self.path.scratch,
                                                "checkpoint.npz")
        # Hidden path to export stats log and figures
        self.path["_optim_output"] = os.path.join(self.path.output, "optimize")

        # Internal check to see if the chosen line search algorithm exists
        if not hasattr(line_search_dir, line_search_method):
            logger.warning(f"{line_search_method} is not a valid line search "
                           f"algorithm, defaulting to 'bracket'")
            line_search_method = "bracket"

        self.line_search_method = line_search_method

        # Internally used parameters for keeping track of optimization
        self._restarted = False
        self._acceptable_preconditioners = ["diagonal"]

        # .title() ensures we grab the class and not the module
        self._line_search = getattr(
            line_search_dir, line_search_method.title())(
                                        step_count_max=step_count_max
                                        )
        
        # Set internal paths to storage locations of models, gradients, etc.
        for tag in ["m_new", "m_old", "m_try",  # Models (vector)
                    "g_new", "g_old", "g_try",  # Gradients (vector)
                    "p_new", "p_old",           # Search Directions (vector)
                    "f_new", "f_old", "f_try",  # Misfits (scalar)
                    "alpha",                    # Step Length (scalar)
                    "dm"                        # Perturbation (vector)
                    "scratch"                   # Scratch array (vector)
                    ]:
            # Make the paths 'hidden' so they don't show up in the parameter 
            # file when we run `seisflows configure`
            self.path[f"_{tag}"] = os.path.join(self.path.scratch, tag)

    def __str__(self):
        """Quickly access underlying line search search history, mostly for
        debug purposes"""
        return self._line_search.__str__()

    @property
    def step_count(self):
        """Convenience property to access `step_count` from line search"""
        return self._line_search.step_count

    def check(self):
        """
        Checks parameters, paths, and dependencies
        """
        if self.preconditioner:
            # This list should match the logic in self.precondition()
            assert self.preconditioner in self._acceptable_preconditioners, \
                f"PRECOND must be in {self._acceptable_preconditioners}"

            assert(os.path.exists(self.path.preconditioner)), (
                f"preconditioner requires PATH.PRECOND pointing to a array-like" 
                f"weight file"
            )

        if self.step_len_init is not None:
            assert 0. < self.step_len_init, \
                f"optimize.step_len_init must be >= 0."
        if self.step_len_max is not None:
            assert 0. < self.step_len_max, \
                f"optimize.step_len_max must be >= 0."
        if self.step_len_min is not None:
            assert 0. < self.step_len_min, \
                f"optimize.step_len_min must be >= 0."
        if self.step_len_max is not None and self.step_len_init is not None:
            assert self.step_len_init < self.step_len_max, \
                f"optimize.step_len_init must be < optimize.step_len_max"
        if self.step_len_min is not None and self.step_len_init is not None:
            assert self.step_len_init > self.step_len_min, \
                f"optimize.step_len_init must be > optimize.step_len_min"
    
        self._line_search.check()

    def setup(self):
        """
        Sets up nonlinear optimization machinery
        """
        unix.mkdir(self.path.scratch)
        unix.mkdir(self.path._optim_output)

        # Load checkpoint (if resuming) or save current checkpoint
        self.load_checkpoint()
        self.checkpoint()  # will be empty

    def checkpoint(self):
        """
        The optimization module (and its underlying `line_search` attribute)
        requires continuity across runs of the same workflow (e.g., in the
        event of a failed job). This function saves internal attributes of
        the optimization module to disk so that a resumed workflow does not
        lose information from its previous version.

        User can checkpoint other variables by adding kwargs
        """
        dict_out = dict(restarted=self._restarted,
                        func_vals=self._line_search.func_vals,
                        step_lens=self._line_search.step_lens,
                        gtg=self._line_search.gtg,
                        gtp=self._line_search.gtp,
                        step_count=self._line_search.step_count,
                        iteridx=self._line_search.iteridx,
                       )
        np.savez(file=self.path._checkpoint, **dict_out)  # NOQA

    def load_checkpoint(self):
        """
        Counterpart to `optimize.checkpoint`. Loads a checkpointed optimization
        module from disk and sets internal tracking attributes.
        """
        # NumPy appends '.npz' when saving. Make sure we honor that.
        if not self.path._checkpoint.endswith(".npz"):
            fid = f"{self.path._checkpoint}.npz"
        else:
            fid = self.path._checkpoint

        if os.path.exists(fid):
            logger.info("re-loading optimization module from checkpoint")
            dict_in = np.load(file=fid)

            self._restarted = bool(dict_in["restarted"])
            self._line_search.func_vals = list(dict_in["func_vals"])
            self._line_search.step_lens = list(dict_in["step_lens"])
            self._line_search.gtg = list(dict_in["gtg"])
            self._line_search.gtp = list(dict_in["gtp"])
            self._line_search.step_count = int(dict_in["step_count"])
            self._line_search.iteridx = list(dict_in["iteridx"])
        else:
            logger.info("no optimization checkpoint file, assume 1st iteration")
            self.checkpoint()

    def precondition(self, q, actions=[], values=[], export_to=None):
        """
        Apply available preconditioner to a given gradient. For efficiency,
        additional `actions` and `values` can be applied while the 
        preconditioner is being applied. This can be called like

        > Gradient.precondition(q=g_new, actions=["*"], values=[-1])

        which would not only apply the preconditioner, but also multiply the
        resulting vector by -1. Multiple actions can be strung together.

        :type q: np.array
        :param q: Vector to precondition, typically gradient contained in: g_new
        :type actions: list of str
        :param actions: additional actions to apply on top of the preconditioner
            such as multiplying by another value
        :type values: list of float
        :param values: values associated with `actions`. Length of values must 
            match actions
        :type export_to: str
        :param export_to: path to export the model values to disk. The filenames
            of the exported files will match the input filenames of this Model,
            not the `other` model. If not given, will export to the same
            directory as the input model `self.path`, which means it overwrites
            the current model
        """
        precon = Model(path=self.path.preconditioner)  # preconditioner

        assert(len(actions) == len(values)), \
            f"number of actions do not match number ofvalues"

        # Apply a diagonal preconditioner, as well as any other user-defined
        # actions ontop of the preconditioner
        if self.preconditioner.upper() == "DIAGONAL":
            logger.info(f"applying diagonal preconditioner and {actions} "
                        f"to {values}")
            q.apply(actions=["*"] + actions, values=[precon] + values,
                    export_to=export_to)

        # > ADD OTHER PRECONDITIONERS HERE, FOLLOW TEMPLATE FOR DIAGONAL
        # if self.preconditioner.upper() == "NAME":
        #   ...

    def compute_direction(self):
        """
        Computes steepest descent search direction (inverse gradient)
        with an optional user-defined preconditioner.

        .. note::
            Other optimization algorithms must overload this method

        :rtype: seisflows.tools.specfem_model.Model
        :return: p_new, search direction as a Model instance
        :raises SystemError: if the search direction is 0, signifying a zero
            gradient which will not do anything in an update
        """
        p_new = Model(path=self.path._g_new)
        if self.preconditioner is not None:
            self.precondition(q=p_new, actions=["*"], values=[-1],
                              export_to=self.path._p_new)
        else:
            # Only multiply by -1 to get the search direction
            p_new.apply(actions=["*"], values=[-1.0], 
                        export_to=self.path._p_new)
            
        # Check the newly created search direction vector 
        p_new = Model(path=self.path._p_new)
        if p_new.get("sum") == 0:
            logger.critical(msg.cli(
                "Search direction vector 'p' is 0, meaning no model update can "
                "take place. Please check your gradient and waveform misfits. "
                "SeisFlows exiting prior to start of line search.", border="=",
                header="optimization gradient error")
            )
            sys.exit(-1)

    def initialize_search(self):
        """
        Setup a the line search machinery by inputting values for the initial
        model. Also contains a check to see if the line search has been 
        restarted.
        """
        f = np.loadtxt(self.path._f_new)  # current misfit value from preprocess
        g = Model(path=self.path._g_new)  # current gradient
        p = Model(path=self.path._p_new)  # current search direction
        
        logger.info("calculating gradient dot products GTG and GTP")
        gtg = g.dot(g)
        gtp = g.dot(p)

        # Restart plugin line search if the optimization library restarts, 
        # restart conditions are determined in `Optimize.compute_direction()`
        if self._restarted:
            self._line_search.clear_search_history()

        # Initialize the line search and save it to disk.
        self._line_search.initialize_line_search(func_val=f, gtg=gtg, gtp=gtp)

    def update_search(self):
        """
        Collect information on a forward evaluation that just took place so 
        that we can assess how to proceed with the current line search. 
        Incremenet the line search step count as we have finished the current.
        """
        # Save step length and associated misfit into the line search
        self._line_search.update_search_history(
            step_len=np.loadtxt(self.path._alpha),  # Step length
            func_val= np.loadtxt(self.path._f_try)  # Misfit for trial model
            )
        logger.info(f"saving misfit and step length for step count == "
                    f"{self.step_count}")

        # Log out the current line search stats for reference
        x, f, idx = self._line_search.get_search_history()
        idx_str = ", ".join([f"{_:>9}" for _ in idx])   # step count
        x_str = ", ".join([f"{_:.3E}" for _ in x])  # step length
        f_str = ", ".join([f"{_:.3E}" for _ in f])  # misfit value
        logger.info(f"step count  = {idx_str}")
        logger.info(f"step length = {x_str}")
        logger.info(f"misfit val  = {f_str}")

        # Increment step count for next line search evaluation
        self._line_search.step_count += 1
        logger.info(f"increment step count -> {self._line_search.step_count}")

    def calculate_step_length(self):
        """
        Determine the step length `alpha` based on the current configuration
        of the line search machinery (i.e., the misfit vs. function evaluation 
        values). Has a catch for whether we want to manual override the step 
        length. Also returns a status that tells the Workflow how to proceed
        with a line search.

        :rtype alpha: float
        :return alpha: step length recommended by the line search. This is
            intrinsically tied to the `status`. When `status`=='TRY', alpha
            represents the step length to take for the next step. When `status`
            =='PASS', then alpha represents the best fitting step count found
            for this line search evaluation.
        :rtype status: str
        :return status: the status recommended to the workflow. Three options:
            Available status returns are:
            'TRY': try/re-try the line search as conditions have not been met
            'PASS': line search was successful, you can terminate the search
            'FAIL': line search has failed for one or more reasons.
        """
        # Used as checks to determine where we are in the inversion
        first_iteration = bool(self._line_search.step_lens.count(0) == 1)
        first_step = bool(self.step_count == 1)
        
        # Initialize Model vectors that will be used for calculating step
        m = Model(path=self.path._m_new)  # current model
        p = Model(path=self.path._p_new)  # current search direction
        norm_m, norm_p = None, None

        # Manually set the step length as some percentage of the model vector.
        # Only do this at the very first model update
        if self.step_len_init and first_iteration and first_step:  
            # We find an initial estimate for step lenght by scaling by the 
            # model divided by the 
            norm_m = m.get(what="absmax")
            norm_p = p.get(what="absmax")
            alpha = self.step_len_init * norm_m / norm_p


            status = None
            logger.debug(f"setting first step length with user-requested "
                         f"`step_len_init`={self.step_len_init}")
        # OR, after the first evaluation, it's expected that the line search 
        # will know how to scale automatically
        else:
            alpha, status = self._line_search.calculate_step_length()

        # Apply optional step length safeguards to prevent step length from 
        # getting too large w.r.t model values
        if status == "TRY" and (self.step_len_max or self.step_len_min):
            # These parameters may have been initiated so we don't do it again
            norm_m = norm_m or m.get(what="absmax")
            norm_p = norm_p or p.get(what="absmax")

            # Determine minimum alpha as a fraction of the current model
            if self.step_len_min:
                logger.debug("checking safeguard min allowable step length: "
                            f"{self.step_len_min * 100}%")
                min_allowable_alpha = self.step_len_min * norm_m / norm_p
                if alpha < min_allowable_alpha:
                    logger.warning(f"step length `alpha` is below minimum "
                                   f"threshold value {min_allowable_alpha:.2f}")
                    if first_step:
                        logger.info("forcing step length to minimum value")
                        alpha = min_allowable_alpha
                    else:
                        logger.critical(msg.mjr("MINIMUM STEP LENGTH EXCEEDED, "
                                                "EXITING WORKFLOW"))
                        sys.exit(-1)
            # Determine maximum alpha as a fraction of the current model
            if self.step_len_max:
                logger.debug("checking safeguard max allowable step length: "
                            f"{self.step_len_max * 100}%")
                max_allowable_alpha = self.step_len_max * norm_m / norm_p
                if alpha > max_allowable_alpha:
                    logger.warning(f"`alpha` has exceeded maximum value "
                                   f"{self.step_len_max * 100}%")
                    if first_step:
                        # Scale back by Golden Ratio so that line search can 
                        # safely increase step length in future steps
                        alpha = 0.618034 * max_allowable_alpha
                        logger.info("reducing step length for first step")
                    else:
                        logger.critical(msg.mjr("MAXIMUM STEP LENGTH EXCEEDED, "
                                                "EXITING WORKFLOW"))
                        sys.exit(-1)

        if alpha is not None:
            logger.info(f"step length `alpha` = {alpha:.3E}")
    
        return alpha, status

    def compute_trial_model(self, alpha):
        """
        Generates a trial model `m_try` by perturbing the starting model `m_new`
        with a given search direction `p_new` and a pre-calculated step length
        `alpha`. Exports `m_try` and `alpha` to disk

        :type alpha: float
        :param alpha: step length recommended by the line search. if None,
            due to failed line search, then this function returns None
        """
        # Current search direction
        p_new = Model(path=self.path._p_new)

        # Get Model Update = step length (a) * step direction (p)
        dm = p_new.apply(actions=["*"], values=[alpha], export_to=self.path._dm)
        logger.info(f"updating model with `dm` "
                    f"(dm_min={dm.get('min'):.2E}, "
                    f"dm_max = {dm.get('max'):.2E})"
                    )
        
        # The new model is the current model `m_new` plus a perturbation `dm`
        # defined as the search direction `p_new` scaled by step length `alpha`
        m_new = Model(path=self.path._m_new)  # This will be overwritten

        # m_i+1 = m_i + a * dm
        m_try = m_new.apply(actions=["+"], values=[dm], 
                            export_to=self.path._m_try)

        # Export `alpha` which is now tied to the current trial model `m_try`
        np.savetxt(self.path._alpha, alpha)

    def finalize_search(self):
        """
        Prepares algorithm machinery and scratch directory for next model update

        - Irretrievably removes `old` model and search parameters
        - Moves current `try` parameters to `old` for reference in next iter.
        - Sets up `new` current parameters
        - Writes and plots statistic output
        - Resets line serach machinery
        """
        unix.cd(self.path.scratch)

        logger.info(msg.sub("FINALIZING LINE SEARCH"))

        # Remove the old-old model parameters (from the last time this was run)
        old_files = glob(os.path.join(self.path.scratch, "?_old*"))
        if old_files:
            logger.info(f"removing previous iteration model files: "
                        f"{old_files}")
            for old_file in old_files:
                unix.rm(old_file)

        # Export stats and figures to output, must run before shifting model
        self.write_stats(fid=os.path.join(self.path._optim_output, 
                                          "output_optim.txt"))
        plot_optim_stats(fid=os.path.join(self.path._optim_output,
                                          "output_optim.txt"),
                         path_out=self.path._optim_output)

        logger.info("setting current model as previous model (?_new -> ?_old)")
        # e.g., f_new.txt -> f_old.txt
        for src in glob(os.path.join(self.path.scratch, "?_new*")):
            # Change 'new' to 'old' in ONLY the filename (#259)
            fid = os.path.basename(src)
            fid_new = fid.replace("_new", "_old")
            dst = src.replace(fid, fid_new)
            unix.mv(src, dst)

        logger.info("setting trial model as starting model (m_try -> m_new)")
        unix.mv(src=self.path._m_try, dst=self.path._m_new)

        # Choose minimum misfit value as final misfit/model, this is the 
        # best accepted model, and does not necessarily correspond to the last
        # line search step
        x, f, idx = self._line_search.get_search_history()

        # Figure out what step length and step count that corresponds to
        np.savetxt(self.path._f_new, f.min())
        logger.info(f"misfit of accepted trial model is f={f.min():.3E}")

        logger.info("resetting line search step count to 0")
        self._line_search.step_count = 0

    def attempt_line_search_restart(self, threshold=1E-3):
        """
        After a failed line search, this determines if restart is worthwhile
        by checking, in effect, if the search direction was the same as the
        negative gradient direction.

        Essentially checking if this is a steepest-descent optimization, which
        cannot and should not be restarted. If the search direction is calc'ed
        by another optimization schema, the search direction and gradient should
        differ and would allow us to restart by falling back to grad. descent.

        :type threshold: float
        :param threshold: angle threshold for the angle between the search
            direction and the gradient.
        :rtype: bool
        :return: pass (True) fail (False) status for retrying line search
        """
        g = Model(path=self.path._g_new)  # Gradient
        p = Model(path=self.path._p_new)  # Search Direction

        # Take the negative gradient and store in scratch vector because we only
        # need it for a short time
        neg_g = g.apply(actions=["*"], values=[-1], 
                        export_to=self.path._scratch)

        # Calculate the angle between the search direction and the neg. gradient
        theta = neg_g.angle(other=p)
        logger.debug(f"checking gradient/search direction angle, "
                     f"theta: {theta:6.3f}")
        
        # Free up scratch
        neg_g.clear()

        if abs(theta) < threshold:
            logger.debug(f"search direction theta={abs(theta):.2E} is below "
                         f"threshold {threshold}, will not attempt restart")
            return False  # Do not restart
        else:
            logger.debug(f"search direction theta={abs(theta):.2E} is above "
                         f"threshold {threshold}, attempting restart")
            return True  # Go for restart

    def restart(self):
        """
        Restarts nonlinear optimization algorithm for any schema that is NOT
        steepest descent (default base class).

        Keeps current position in model space, but discards history of
        nonlinear optimization algorithm in an attempt to recover from
        numerical stagnation.

        .. note::
            steepest descent optimization algorithm does not have any restart
            capabilities. This function is instantiated here to be overwritten
            by child classes
        """
        pass

    def write_stats(self, fid="./optim_optim.txt"):
        """
        Write stats to file so that we don't lose information to subsequent
        iterations. File is written to path._output_optim/fid

        :type fid: str
        :param fid: full path and filename to save the text file. defaults to
            ./output_optim.txt
        """
        logger.info(f"writing optimization stats: '{fid}'")

        keys = ["misfit", "step_count",  "step_length", 
                "slope", "theta", "if_restarted"
                ]

        # First time, write header information and start model misfit. Note that
        # most of the statistics do not apply to the starting model so they
        # are set to 0 by default
        if not os.path.exists(fid):
            with open(fid, "w") as f_:
                x, f, _  = self._line_search.get_search_history()
                header = ",".join(keys) + "\n"
                f_.write(header)
                # Write values for first iteration
                _write_vals = []
                for key in keys:
                    if key == "step_length":
                        val = x[0]
                    elif key == "misfit":
                        val = f[0]
                    else:
                        val = 0
                    _write_vals.append(f"{val:6.3E}") 
                stats_str = ",".join(_write_vals) + "\n"
                f_.write(stats_str)  

        # Write stats for the current, finished, line search
        stats = self._get_stats()
        with open(fid, "a") as f_:
            stats_str = [f"{stats[key]:6.3E}" for key in keys]
            stats_str = ",".join(stats_str) + "\n"
            f_.write(stats_str)

    def _get_stats(self):
        """
        Get Optimization statistics for the current evaluation of an inversion.
        Returns a dictionary of values which can then be written to a text file
        or plotted for convenience.           

        .. note::

            The following `factor` was included in the original SeisFlows stats
            but started throwing runtime errors. Not sure what is was for -- BC

            factor = -1 * dot(g.vector, g.vector)
            factor = factor ** -0.5 * (f[1] - f[0]) / (x[1] - x[0]) 

        .. note:: L1/L2 Norm

            BC The original SeisFlows calculated L1 and L2 norm of the gradient.
            I removed these because they are heavy operations to run on large
            model vectors and were not very useful when evaluating models. 
            Maybe okay for 2D problems where you can quickly take a model norm
            but large 3D models were hanging a long time. If we want to 
            re-instate them then we need to find a more efficient way to 
            calculate norms with the new Model class paradigm (#245)
                    
        :rtype: Dict of float
        :return: SeisFlows dictionary object containing statistical informaion
            about the optimization module
        """
        # We construct stats from the gradient (g), search direction (p),
        # step lengths (x) and function values/misfit (f)
        g = Model(self.path._g_new)
        p = Model(self.path._p_new)

        x, f, _  = self._line_search.get_search_history()
        step_count = self._line_search.step_count

        # Construct statistics
        step_length = x[f.argmin()]  # final, accepted step length
        misfit = f[f.argmin()]  # final, accepted misfit value
        if_restarted = int(self._restarted)
        slope = (f[1] - f[0]) / (x[1] - x[0])  # Slope of the misfit function

        # Calculate norms which may take a little bit of time
        # !!! See note above
        # grad_norm_L1 = np.linalg.norm(g.vector, 1)  # L1 norm of gradient
        # grad_norm_L2 = np.linalg.norm(g.vector, 2)  # L2 norm of gradient

        # Deviation of the search direction from the gradient direction
        # Take the negative gradient and store in scratch vector
        neg_g = g.apply(actions=["*"], values=[-1], 
                        export_to=self.path._scratch)
        theta = 180. * np.pi ** -1 * p.angle(neg_g) 

        dict_out = Dict(step_count=step_count, step_length=step_length,
                        misfit=misfit, if_restarted=if_restarted,
                        slope=slope, theta=theta)
        
        # Clear out disk space
        neg_g.clear()
    
        return dict_out
