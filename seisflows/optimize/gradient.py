#!/usr/bin/env python3
"""
Gradient descent nonlinear optimization algorithm. Acts as the Base class for
optimization.

The Optimization library contains classes and methods used to solve nonlinear
optimization problems, i.e., misfit minimization. Various subclasses implement
different optimization algorithms.

.. note::
    By default the base class implements a steepest descent optimization
"""
import os
import sys
import logging
import numpy as np

from seisflows.tools import msg, unix
from seisflows.tools.math import angle, dot
from seisflows.plugins import line_search, preconds
from seisflows.tools.math import poissons_ratio
from seisflows.config import SeisFlowsPathsParameters

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]
solver = sys.modules["seisflows_solver"]


class Gradient:
    """
    Nonlinear optimization abstract base class poviding a gradient/steepest
    descent optimization algorithm.

    Nonlinear conjugate, quasi-Newton and Newton methods can be implemented on
    top of this base class.

    .. note::
        To reduce memory overhead, vectors are read from disk rather than passed
        from calling routines. For example, at the beginning of
        compute_direction the current gradient is read from  'g_new' and the
        resulting search direction is written to 'p_new'. As the inversion
        progresses, other information is stored as well.

    .. note::
        The default numerical parameters defined below should work well for a
        range of applications without manual tuning. If the nonlinear
        optimization procedure stagnates, it may be due to issues involving
        data quality or the choice of data misfit, data processing, or
        regularization parameters.  Problems in any of these areas usually
        manifest themselves through stagnation of the nonlinear optimization
        algorithm.

    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        Initialize internally used variables for optimization workflow

        :type iter: int
        :param iter: the current iteration of the workflow
        :type line_search: Class
        :param line_search: a class controlling the line search functionality
            for determining step length
        :type precond: Class
        :param precond: a class controlling the preconditioner functionality
            for preconditiong gradient information
        :type restarted: bool
        :param restarted: a flag signalling if the optimization algorithm has
            been restarted recently
        """
        self.iter = 1
        self.line_search = None
        self.precond = None
        self.restarted = False

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        # Define the Parameters required by this module
        sf.par("LINESEARCH", required=False, default="Bracket", par_type=str,
               docstr="Algorithm to use for line search, see "
                      "seisflows.plugins.line_search for available choices")

        sf.par("PRECOND", required=False, par_type=str,
               docstr="Algorithm to use for preconditioning gradients, see "
                      "seisflows.plugins.preconds for available choices")

        sf.par("STEPCOUNTMAX", required=False, default=10, par_type=int,
               docstr="Max number of trial steps in line search before a "
                      "change in line search behavior")

        sf.par("STEPLENINIT", required=False, default=0.05, par_type=float,
               docstr="Initial line search step length, as a fraction "
                      "of current model parameters")

        sf.par("STEPLENMAX", required=False, default=0.5, par_type=float,
               docstr="Max allowable step length, as a fraction of "
                      "current model parameters")

        # Define the Paths required by this module
        sf.path("OPTIMIZE", required=False,
                default=os.path.join(PATH.SCRATCH, "optimize"),
                docstr="scratch path for nonlinear optimization data")

        return sf

    def check(self, validate=True):
        """
        Checks parameters, paths, and dependencies
        """
        if validate:
            self.required.validate()

        if PAR.LINESEARCH:
            assert PAR.LINESEARCH in dir(line_search), \
                f"LINESEARCH parameter must be in {dir(line_search)}"

        if PAR.PRECOND:
            assert PAR.PRECOND in dir(preconds), \
                f"PRECOND must be in {dir(preconds)}"

        assert 0. < PAR.STEPLENINIT, f"STEPLENINIT must be >= 0."
        assert 0. < PAR.STEPLENMAX, f"STEPLENMAX must be >= 0."
        assert PAR.STEPLENINIT < PAR.STEPLENMAX, \
            f"STEPLENINIT must be < STEPLENMAX"

    def setup(self, m_new=None):
        """
        Sets up nonlinear optimization machinery

        :type m_new: np.array
        :param m_new: Initial model vector which can be user-provided. Otherwise
            the initial model will be read from PATH.MODEL_INIT
        """
        unix.mkdir(PATH.OPTIMIZE)

        # Line search machinery is defined externally as a plugin class
        if PAR.LINESEARCH:
            self.line_search = getattr(line_search, PAR.LINESEARCH)(
                step_count_max=PAR.STEPCOUNTMAX,
                step_len_max=PAR.STEPLENMAX,
            )
        if PAR.PRECOND:
            self.precond = getattr(preconds, PAR.PRECOND)()

        # Read in initial model as a vector and ensure it is a valid model
        if m_new is None:
            assert(os.path.exists(PATH.MODEL_INIT)), (
                "optimization library requires that PATH.MODEL_INIT exists"
            )
            m_new = solver.merge(solver.load(PATH.MODEL_INIT))

        np.save(self.vectors("m_new"), m_new)
        self.check_model(m_new)

    @property
    def eval_str(self):
        """
        Print out the evaluation string, which states what iteration and line
        search step count we are at. Useful for log statements

        For example, an inversion at iteration 1 and step count 2 will return
        'i01s02'
        """
        iter_ = self.iter
        step = self.line_search.step_count
        return f"i{iter_:0>2}s{step:0>2}"

    def vectors(self, name):
        """
        Convenience function to access the full paths of model and gradient
        vectors that are saved to disk

        .. note:: the available options that can be created
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
        # Set model and gradient filenames as attributes
        acceptable_names = ["m_new", "m_old", "m_try",
                            "g_new", "g_old", "g_try",
                            "p_new", "p_old", "alpha",
                            "f_new", "f_old", "f_try"]
        assert(name in acceptable_names)
        return os.path.join(PATH.OPTIMIZE, f"{name}.npy")

    def compute_direction(self):
        """
        Computes a steepest descent search direction (inverse gradient)
        with an optional user-defined preconditioner.

        .. note::
            Other optimization algorithms must overload this method
        """
        self.logger.info(f"computing search direction with {PAR.OPTIMIZE}")

        g_new = np.load(self.vectors("g_new"))
        if self.precond is not None:
            p_new = -1 * self.precond(g_new)
        else:
            p_new = -1 * g_new
        np.save(self.vectors("p_new"), p_new)

    def initialize_search(self):
        """
        Initialize the plugin line search machinery. Should only be run at
        the beginning of line search, by the main workflow module.
        """
        m = np.load(self.vectors("m_new"))
        g = np.load(self.vectors("g_new"))
        p = np.load(self.vectors("p_new"))
        f = np.load(self.vectors("f_new"))

        norm_m = max(abs(m))
        norm_p = max(abs(p))

        gtg = dot(g, g)
        gtp = dot(g, p)

        # Restart plugin line search if the optimization library restarts
        if self.restarted:
            self.line_search.clear_history()

        # Optional safeguard to prevent step length from getting too large
        if PAR.STEPLENMAX:
            self.line_search.step_len_max = PAR.STEPLENMAX * norm_m / norm_p
            self.logger.debug(f"max step length safeguard is: "
                              f"{self.line_search.step_len_max:.2E}")

        self.line_search.initialize(step_len=0., func_val=f, gtg=gtg, gtp=gtp)
        alpha, _ = self.line_search.calculate_step()

        # Alpha defines the trial step length. Optional step length override
        if PAR.STEPLENINIT and len(self.line_search.step_lens) <= 1:
            alpha = PAR.STEPLENINIT * norm_m / norm_p
            self.logger.debug(f"overwrite initial step length: {alpha:.2E}")

        # The new model is the old model, scaled by the step direction and
        # gradient threshold to remove any outlier values
        m_try = m + alpha * p

        np.save(self.vectors("m_try"), m_try)
        np.save(self.vectors("alpha"), alpha)
        self.check_model(m_try)

    def update_search(self):
        """
        Updates line search status and step length and checks if the line search
        has been completed.

        Available status codes from line_search.update():
            status == 1  : finished
            status == 0 : not finished
            status == -1  : failed
        """
        self.line_search.update(step_len=np.load(self.vectors("alpha")),
                                func_val=np.load(self.vectors("f_try")))
        alpha, status = self.line_search.calculate_step()

        # New search direction needs to be searchable on disk
        if status in [0, 1]:
            m = np.load(self.vectors("m_new"))
            p = np.load(self.vectors("p_new"))
            np.save(self.vectors("alpha"), alpha)

            m_try = m + alpha * p
            np.save(self.vectors("m_try"), m_try)
            self.check_model(m_try)

        return status

    def finalize_search(self):
        """
        Prepares algorithm machinery and scratch directory for next model update

        Removes old model/search parameters, moves current parameters to old,
        sets up new current parameters and writes statistic outputs
        """
        unix.cd(PATH.OPTIMIZE)
        self.logger.info(msg.sub("FINALIZING LINE SEARCH"))

        # Remove the old model parameters
        if self.iter > 1:
            self.logger.info("removing previously accepted model files (old)")
            for fid in [self.vectors("m_old"), self.vectors("f_old"),
                        self.vectors("g_old"), self.vectors("p_old")]:
                unix.rm(fid)

        # Needs to be run before shifting model in next step
        self.write_stats()

        self.logger.info("shifting current model (new) to previous model (old)")
        unix.mv(self.vectors("m_new"), self.vectors("m_old"))
        unix.mv(self.vectors("f_new"), self.vectors("f_old"))
        unix.mv(self.vectors("g_new"), self.vectors("g_old"))
        unix.mv(self.vectors("p_new"), self.vectors("p_old"))

        self.logger.info("setting accepted line search model as current model")
        unix.mv(self.vectors("m_try"), self.vectors("m_new"))

        f = self.line_search.search_history()[1]
        np.save(self.vectors("f_new"), f.min())
        self.logger.info(f"current misfit is {f.min():.3E}")

        self.logger.info("resetting line search step count to 0")
        self.line_search.step_count = 0

    def retry_status(self):
        """
        After a failed line search, this determines if restart is worthwhile
        by checking, in effect, if the search direction was the same as gradient
        direction
        """
        g = np.load(self.vectors("g_new"))
        p = np.load(self.vectors("p_new"))
        theta = angle(p, -g)

        self.logger.debug(f"theta: {theta:6.3f}")

        thresh = 1.e-3
        if abs(theta) < thresh:
            return 0
        else:
            return 1

    def restart(self):
        """
        Restarts nonlinear optimization algorithm for any schema that is NOT
        steepest descent (default base class).

        Keeps current position in model space, but discards history of
        nonlinear optimization algorithm in an attempt to recover from
        numerical stagnation.
        """
        # Steepest descent (base) does not need to be restarted
        if PAR.OPTIMIZE.capitalize() != "Gradient":
            g = np.load(self.vectors("g_new"))
            np.save(self.vectors("p_new"), -g)

            self.line_search.clear_history()
            self.restarted = 1

    def write_stats(self):
        """
        Simplified write function to append values to text files.
        Used because stats line search information can be overwritten
        by subsequent iterations so we need to append values to text files
        if they should be retained.

        :type log: str
        :param log: name of the file to write to. Will append .txt to it
        :type value: float
        :param value: value to write to file
        :type format: str
        :param format: string formatter for value
        """
        self.logger.info(f"writing optimization stats")
        fid = os.path.join(PATH.OUTPUT,  f"optim_stats.txt")

        # First time, write header information
        if not os.path.exists(fid):
            with open(fid, "w") as f:
                for header in ["ITER", "FACTOR", "GRAD_NORM_L1", "GRAD_NORM_L2",
                              "MISFIT", "RESTART", "SLOPE", "STEP", "LENGTH",
                              "THETA"]:
                    f.write(f"{header.upper()},")
                f.write("\n")

        g = np.load(self.vectors("g_new"))
        p = np.load(self.vectors("p_new"))
        x = self.line_search.search_history()[0]
        f = self.line_search.search_history()[1]

        # Calculated stats factors
        factor = -dot(g, g) ** -0.5 * (f[1] - f[0]) / (x[1] - x[0])
        grad_norm_L1 = np.linalg.norm(g, 1)
        grad_norm_L2 = np.linalg.norm(g, 2)
        misfit = f[0]
        restarted = self.restarted
        slope = (f[1] - f[0]) / (x[1] - x[0])
        step_count = self.line_search.step_count
        step_length = x[f.argmin()]
        theta = 180. * np.pi ** -1 * angle(p, -g)

        with open(fid, "a") as f:
            f.write(f"{self.iter:0>2},"
                    f"{factor:6.3E},"
                    f"{grad_norm_L1:6.3E},"
                    f"{grad_norm_L2:6.3E},"
                    f"{misfit:6.3E},"
                    f"{restarted:6.3E},"
                    f"{slope:6.3E},"
                    f"{step_count:0>2},"
                    f"{step_length:6.3E},"
                    f"{theta:6.3E}\n"
                    )

    def check_model(self, m, min_pr=-1., max_pr=0.5):
        """
        Check to ensure that the model parameters fall within the guidelines
        of the solver. Print off min/max model parameters for the User.

        :type m: np.array
        :param m: model to check parameters of
        :type min_pr: float
        :param min_pr: minimum allowable Poisson's ratio value dictated by
            SPECFEM
        :type max_pr: float
        :param max_pr: maximum allowable Poisson's ratio value dictated by
            SPECFEM
        """
        # Dynamic way to split up the model based on number of params
        pars = {}
        for i, par in enumerate(solver.parameters):
            pars[par] = np.split(m, len(solver.parameters))[i]

        # Check Poisson's ratio, which will error our SPECFEM if outside limits
        if (pars["vp"] is not None) and (pars["vs"] is not None):
            self.logger.debug(f"checking poissons ratio")
            pars["pr"] = poissons_ratio(vp=pars["vp"], vs=pars["vs"])
            if pars["pr"].min() < 0:
                self.logger.warning("minimum poisson's ratio is negative")
            if pars["pr"].min() < min_pr:
                self.logger.warning(f"minimum poisson's ratio out of bounds: "
                                    f"{pars['pr'].max()} > {max_pr}")
            if pars["pr"].max() > max_pr:
                self.logger.warning(f"maximum poisson's ratio out of bounds: " 
                                    f"{pars['pr'].min()} < {min_pr}")

        # Tell the User min and max values of the updated model
        self.logger.info(f"model parameters")
        parts = "{minval:.2f} <= {key} <= {maxval:.2f}"
        for key, vals in pars.items():
            self.logger.info(parts.format(minval=vals.min(), key=key,
                                          maxval=vals.max())
                             )






