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
from seisflows.tools.math import angle, dot
from seisflows.tools.model import Model
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
    :param step_len_max: maximum allowable step length during the line
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
                 preconditioner=None, step_count_max=10, step_len_init=0.05,
                 step_len_max=0.5, workdir=os.getcwd(), path_optimize=None,
                 path_output=None, path_preconditioner=None, **kwargs):
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
        self._acceptable_vectors = ["m_new", "m_old", "m_try",
                                    "g_new", "g_old", "g_try",
                                    "p_new", "p_old", "alpha",
                                    "f_new", "f_old", "f_try"]
        self._acceptable_preconditioners = ["diagonal"]

        # .title() ensures we grab the class and not the module
        self._line_search = getattr(
            line_search_dir, line_search_method.title())(
            step_count_max=step_count_max, step_len_max=step_len_max,
        )

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
        if self.step_len_max is not None and self.step_len_init is not None:
            assert self.step_len_init < self.step_len_max, \
                f"optimize.step_len_init must be < optimize.step_len_max"
    
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

    def load_vector(self, name):
        """
        Convenience function to access the full paths of model and gradient
        vectors that are saved to disk

        .. note::
            Model, gradient and misfit vectors are named as follows:
            m_new: current model
            m_old: previous model
            m_try: line search model
            f_new: current objective function value
            f_old: previous objective function value
            f_try: line search function value
            g_new: current gradient direction
            g_old: previous gradient direction
            p_new: current search direction
            p_old: previous search direction
            alpha: trial search direction (aka p_try)

        :type name: str
        :param name: name of the vector, acceptable: m, g, p, f, alpha
        """
        assert(name in self._acceptable_vectors)

        model_npz = os.path.join(self.path.scratch, f"{name}.npz")
        model_npy = model_npz.replace(".npz", ".npy")
        model_txt = model_npz.replace(".npz", ".txt")

        if os.path.exists(model_npz):
            model = Model(path=model_npz)
        elif os.path.exists(model_npy):
            model = np.load(model_npy)
        elif os.path.exists(model_txt):
            model = float(np.loadtxt(model_txt))
        else:
            raise FileNotFoundError(f"no optimization file found for '{name}'")

        return model

    def save_vector(self, name, m):
        """
        Convenience function to save/overwrite vectors on disk

        :type name: str
        :param name: name of the vector to overwrite
        :type m: seisflows.tools.specfem.Model or float
        :param m: Model to save to disk as npz array
        """
        assert(name in self._acceptable_vectors)

        if isinstance(m, Model):
            path = os.path.join(self.path.scratch, f"{name}.npz")
            m.model = m.split()  # overwrite m representation
            m.save(path=path)
        elif isinstance(m, np.ndarray):
            path = os.path.join(self.path.scratch, f"{name}.npy")
            np.save(path=path)
        elif isinstance(m, (float, int)):
            path = os.path.join(self.path.scratch, f"{name}.txt")
            np.savetxt(path, [m])
        else:
            raise TypeError(f"optimize.save unrecognized type error {type(m)}")

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
                        step_len_max=self._line_search.step_len_max,
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
            self._line_search.step_len_max = float(dict_in["step_len_max"])
            self._line_search.iteridx = list(dict_in["iteridx"])
        else:
            logger.info("no optimization checkpoint file, assume 1st iteration")
            self.checkpoint()

    def _precondition(self, q):
        """
        Apply available preconditioner to a given gradient

        :type q: np.array
        :param q: Vector to precondition, typically gradient contained in: g_new
        :rtype: np.array
        :return: preconditioned vector
        """
        if self.preconditioner is not None:
            p = Model(path=self.path.preconditioner)
            if self.preconditioner.upper() == "DIAGONAL":
                logger.info("applying diagonal preconditioner")
                return p.vector * q
            else:
                raise NotImplementedError(
                    f"preconditioner {self.preconditioner} not supported"
                )
        else:
            return q

    def compute_direction(self):
        """
        Computes steepest descent search direction (inverse gradient)
        with an optional user-defined preconditioner.

        .. note::
            Other optimization algorithms must overload this method

        :rtype: seisflows.tools.specfem.Model
        :return: search direction as a Model instance
        :raises SystemError: if the search direction is 0, signifying a zero
            gradient which will not do anything in an update
        """
        g_new = self.load_vector("g_new")
        p_new = g_new.copy()
        p_new.update(vector=-1 * self._precondition(g_new.vector))

        if sum(p_new.vector) == 0:
            logger.critical(msg.cli(
                "Search direction vector 'p' is 0, meaning no model update can "
                "take place. Please check your gradient and waveform misfits. "
                "SeisFlows exiting prior to start of line search.", border="=",
                header="optimization gradient error")
            )
            sys.exit(-1)

        return p_new

    def initialize_search(self):
        """
        Setup a the line search machinery by inputting values for the initial
        model. Also contains a check to see if the line search has been 
        restarted.
        """
        f = self.load_vector("f_new")  # current misfit value from preprocess
        g = self.load_vector("g_new")  # current gradient from scaled kernels
        p = self.load_vector("p_new")  # current search direction
        gtg = dot(g.vector, g.vector)
        gtp = dot(g.vector, p.vector)

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
        alpha_try = self.load_vector("alpha")  # step length
        f_try = self.load_vector("f_try")  # misfit for the trial model
        self._line_search.update_search_history(step_len=alpha_try, 
                                                func_val=f_try)
        logger.info(f"saving misfit and step length for step count == "
                    f"{self.step_count}")

        # Log out the current line search stats for reference
        x, f = self._line_search.get_search_history()
        i_str = ", ".join([f"{_:>9}" for _ in range(len(x))])
        x_str = ", ".join([f"{_:.3E}" for _ in x])
        f_str = ", ".join([f"{_:.3E}" for _ in f])
        logger.info(f"step count  = {i_str}")
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
        
        # Initialize empty vectors for checks
        m, p = None, None

        # Manually set the step length as some percentage of the model vector.
        # Only do this at the very first model update
        if self.step_len_init and first_iteration and first_step:
            m = self.load_vector("m_new")  # current model
            p = self.load_vector("p_new")  # current search direction
            norm_m = max(abs(m.vector))
            norm_p = max(abs(p.vector))

            alpha = self.step_len_init * norm_m / norm_p
            status = None
            logger.debug(f"setting first step length with user-requested "
                         f"`step_len_init`={self.step_len_init}")
        # OR, after the first evaluation, it's expected that the line search 
        # will know how to scale automatically
        else:
            alpha, status = self._line_search.calculate_step_length()

        # OPTIONAL: Apply step length safeguard to prevent step length from 
        # getting too large w.r.t model values
        if status == "TRY" and self.step_len_max:
            logger.debug("checking safeguard: maximum allowable step length")
            m = m or self.load_vector("m_new")  # current model
            p = p or self.load_vector("p_new")  # current search direction
            norm_m = max(abs(m.vector))
            norm_p = max(abs(p.vector))

            # Determine maximum alpha as a fraction of the current model
            max_allowable_alpha = self.step_len_max * norm_m / norm_p
            if alpha > max_allowable_alpha:
                logger.warning(f"safeguard: alpha has exceeded maximum step "
                               f"length {self.step_len_max}, capping value")
                if first_step:
                    # If this is the first step, pull back slightly so that line 
                    # search can safely increase step length later
                    alpha = 0.618034 * max_allowable_alpha
                else:
                    alpha = max_allowable_alpha

        if alpha is not None:
            logger.info(f"step length `alpha` = {alpha:.3E}")
    
        return alpha, status

    def compute_trial_model(self, alpha):
        """
        Generates a trial model `m_try` by perturbing the starting model `m_new`
        with a given search direction `p_new` and a pre-calculated step length
        `alpha`

        :type alpha: float
        :param alpha: step length recommended by the line search. if None,
            due to failed line search, then this function returns None
        :rtype: np.array or NoneType
        :return: trial model that can be used for line search evaluation or 
            None if alpha is None
        """
        if alpha is not None:
            # The new model is the old model plus a step with a given magnitude 
            m_try = self.load_vector("m_new").copy()
            p = self.load_vector("p_new")  # current search direction

            dm = alpha * p.vector  # update = step length * step direction
            logger.info(f"updating model with `dm` (dm_min={dm.min():.2E}, "
                        f"dm_max = {dm.max():.2E})")
            m_try.update(vector=m_try.vector + dm)
        else:
            m_try = None

        return m_try

    def finalize_search(self):
        """
        Prepares algorithm machinery and scratch directory for next model update

        Removes old model/search parameters, moves current parameters to old,
        sets up new current parameters and writes statistic outputs
        """
        unix.cd(self.path.scratch)

        logger.info(msg.sub("FINALIZING LINE SEARCH"))

        # Remove the old-old model parameters (from the last time this was run)
        if glob("?_old*"):
            logger.info("removing previously accepted model files (?_old)")
            for fid in ["m_old.npz", "f_old.txt", "g_old.npz", "p_old.npz"]:
                unix.rm(os.path.join(self.path.scratch, fid))

        # Export stats and figures to output, must run before shifting model
        self.write_stats(fid=os.path.join(self.path._optim_output, 
                                          "output_optim.txt"))
        plot_optim_stats(fid=os.path.join(self.path._optim_output,
                                          "output_optim.txt"),
                         path_out=self.path._optim_output)

        logger.info("setting current model as previous model (new -> old)")
        # e.g., m_new.npz -> m_old.npz
        for src in glob(os.path.join(self.path.scratch, "*_new.*")):
            dst = src.replace("_new.", "_old.")
            unix.mv(src, dst)

        logger.info("setting trial model as starting model (m_try -> m_new)")
        unix.mv(src=os.path.join(self.path.scratch, "m_try.npz"),
                dst=os.path.join(self.path.scratch, "m_new.npz"))

        # Choose minimum misfit value as final misfit/model. index 0 is initial
        x, f = self._line_search.get_search_history()
        self.save_vector("f_new", f.min())
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
        differ

        :type threshold: float
        :param threshold: angle threshold for the angle between the search
            direction and the gradient.
        :rtype: bool
        :return: pass (True) fail (False) status for retrying line search
        """
        g = self.load_vector("g_new")
        p = self.load_vector("p_new")

        theta = angle(p.vector, -1 * g.vector)
        logger.debug(f"checking gradient/search direction angle, "
                     f"theta: {theta:6.3f}")

        if abs(theta) < threshold:
            logger.info(f"search direction below threshold {threshold}, will "
                        f"not attempt restart")
            return False  # Do not restart
        else:
            logger.info(f"search direction above threshold {threshold}, "
                        f"attempting restart")
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

    def get_stats(self):
        """
        Get Optimization statistics for the current evaluation of an inversion.
        Returns a dictionary of values which can then be written to a text file
        or plotted for convenience.           

        .. note::

            The following `factor` was included in the original SeisFlows stats
            but started throwing runtime errors. Not sure what is was for -- BC

            factor = -1 * dot(g.vector, g.vector)
            factor = factor ** -0.5 * (f[1] - f[0]) / (x[1] - x[0]) 
                    
        :rtype: Dict of float
        :return: SeisFlows dictionary object containing statistical informaion
            about the optimization module
        """
        # We construct stats from the gradient (g), search direction (p),
        # step lengths (x) and function values/misfit (f)
        g = self.load_vector("g_new")
        p = self.load_vector("p_new")
        x, f  = self._line_search.get_search_history()
        step_count = self._line_search.step_count

        # Construct statistics
        step_length = x[f.argmin()]  # final, accepted step length
        misfit = f[f.argmin()]  # final, accepted misfit value
        if_restarted = int(self._restarted)
        grad_norm_L1 = np.linalg.norm(g.vector, 1)  # L1 norm of gradient
        grad_norm_L2 = np.linalg.norm(g.vector, 2)  # L2 norm of gradient
        slope = (f[1] - f[0]) / (x[1] - x[0])  # Slope of the misfit function
        # Deviation of the search direction from the gradient direction
        theta = 180. * np.pi ** -1 * angle(p.vector, -1 * g.vector)  

        dict_out = Dict(step_count=step_count, step_length=step_length,
                        misfit=misfit, if_restarted=if_restarted,
                        grad_norm_L1=grad_norm_L1, grad_norm_L2=grad_norm_L2, 
                        slope=slope, theta=theta)
    
        return dict_out

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
                "grad_norm_L1", "grad_norm_L2",
                "slope", "theta", "if_restarted"
                ]

        # First time, write header information and start model misfit. Note that
        # most of the statistics do not apply to the starting model so they
        # are set to 0 by default
        if not os.path.exists(fid):
            with open(fid, "w") as f_:
                x, f  = self._line_search.get_search_history()
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
        stats = self.get_stats()
        with open(fid, "a") as f_:
            stats_str = [f"{stats[key]:6.3E}" for key in keys]
            stats_str = ",".join(stats_str) + "\n"
            f_.write(stats_str)

