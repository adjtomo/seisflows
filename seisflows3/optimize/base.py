#!/usr/bin/env python3
"""
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

from seisflows3.tools import msg, unix
from seisflows3.tools.math import angle, dot
from seisflows3.plugins import line_search, preconds
from seisflows3.tools.specfem import check_poissons_ratio
from seisflows3.config import SeisFlowsPathsParameters, CFGPATHS

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]
solver = sys.modules["seisflows_solver"]


class Base:
    """
    Nonlinear optimization abstract base class.

    This base class provides a steepest descent optimization algorithm.

    Nonlinear conjugate, quasi-Newton and Newton methods can be implemented on
    top of this base class.

    .. note::
        To reduce memory overhead, vectors are read from disk rather than passed
        from calling routines. For example, at the beginning of
        compute_direction the current gradient is  read from  'g_new' and the
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
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings

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

        :param m_new: current model
        :param m_old: previous model
        :param m_try: line search model
        :param f_new: current objective function value
        :param f_old: previous objective function value
        :param f_try: line search function value
        :param g_new: current gradient direction
        :param g_old: previous gradient direction
        :param p_new: current search direction
        :param p_old: previous search direction
        """
        self.iter = 1
        self.line_search = None
        self.precond = None
        self.restarted = False

        # Define the names of output stats logs to keep all paths in one place
        # Line search log is named differently so that optimize doesn't
        # overwrite this log file when intiating the stats directory
        self.line_search_log = "line_search"
        self.log_factor = "factor"
        self.log_gradient_norm_L1 = "gradient_norm_L1"
        self.log_gradient_norm_L2 = "gradient_norm_L2"
        self.log_misfit = "misfit"
        self.log_restarted = "restarted"
        self.log_slope = "slope"
        self.log_step_count = "step_count"
        self.log_step_length = "step_length"
        self.log_theta = "theta"

        # Define the names of variables used to keep track of models etc. so
        # that we don't have multiple strings floating around defining the same
        # thing
        self.m_new = "m_new.npy"
        self.m_old = "m_old.npy"
        self.m_try = "m_try.npy"
        self.f_new = "f_new.npy"
        self.f_old = "f_old.npy"
        self.f_try = "f_try.npy"
        self.g_new = "g_new.npy"
        self.g_old = "g_old.npy"
        self.p_new = "p_new.npy"
        self.p_old = "p_old.npy"
        self.alpha = "alpha.npy"

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
                      "seisflows3.plugins.line_search for available choices")

        sf.par("PRECOND", required=False, par_type=str,
               docstr="Algorithm to use for preconditioning gradients, see "
                      "seisflows3.plugins.preconds for available choices")

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
        self.logger.debug(msg.check(type(self)))

        if validate:
            self.required.validate()

        if PAR.OPTIMIZE == "base":
            print(msg.cli("Base optimization cannot be used standalone, and "
                          "must be over-loaded by a subclass. You can run"
                          "'seisflows print module' to find available "
                          "subclasses", header="error", border="="))
            sys.exit(-1)

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

    def setup(self):
        """
        Sets up nonlinear optimization machinery
        """
        msg.setup(type(self))

        # All ptimization statistics text files will be written to path_stats
        path_stats = os.path.join(PATH.WORKDIR, CFGPATHS.STATSDIR)
        unix.mkdir(path_stats)

        # Line search machinery is defined externally as a plugin class
        self.line_search = getattr(line_search, PAR.LINESEARCH)(
            step_count_max=PAR.STEPCOUNTMAX, step_len_max=PAR.STEPLENMAX,
            log_file=os.path.join(path_stats, f"{self.line_search_log}.txt"),
        )

        if PAR.PRECOND:
            self.precond = getattr(preconds, PAR.PRECOND)()
        else:
            self.precond = None

        # Instantiate all log files in stats/ directory as empty text files
        # OVERWRITES any existing stats/ log files that may already be there
        for key, val in vars(self).items():
            if "log_" in key:
                self.write_stats(val)

        # Ensure that line search step count starts at 0 (workflow.intialize)
        self.write_stats(self.log_step_count, 0)

        unix.mkdir(PATH.OPTIMIZE)
        if "MODEL_INIT" in PATH:
            m_new = solver.merge(solver.load(PATH.MODEL_INIT))
            self.save(self.m_new, m_new)
            self.check_model(m_new, self.m_new)

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

    def compute_direction(self):
        """
        Computes a steepest descent search direction (inverse gradient)
        with an optional user-defined preconditioner.

        .. note::
            Other optimization algorithms must overload this method
        """
        self.logger.info(f"computing search direction with {PAR.OPTIMIZE}")

        g_new = self.load(self.g_new)
        if self.precond is not None:
            p_new = -1 * self.precond(g_new)
        else:
            p_new = -1 * g_new
        self.save(self.p_new, p_new)

    def initialize_search(self):
        """
        Initialize the plugin line search machinery. Should only be run at
        the beginning of line search, by the main workflow module.
        """
        m = self.load(self.m_new)
        g = self.load(self.g_new)
        p = self.load(self.p_new)
        f = self.loadtxt(self.f_new)
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

        # Alpha defines the trial step length
        alpha, _ = self.line_search.initialize(iter=self.iter, step_len=0.,
                                               func_val=f, gtg=gtg, gtp=gtp
                                               )

        # Optional initial step length override
        if PAR.STEPLENINIT and len(self.line_search.step_lens) <= 1:
            alpha = PAR.STEPLENINIT * norm_m / norm_p
            self.logger.debug(f"manually set initial step length: {alpha:.2E}")

        # The new model is the old model, scaled by the step direction and
        # gradient threshold to remove any outlier values
        m_try = m + alpha * p

        self.save(self.m_try, m_try)
        self.savetxt(self.alpha, alpha)
        self.check_model(m_try, self.m_try)

    def update_search(self):
        """
        Updates line search status and step length and checks if the line search
        has been completed.

        Available status codes from line_search.update():
            status == 1  : finished
            status == 0 : not finished
            status == -1  : failed
        """
        alpha, status = self.line_search.update(
            iter=self.iter, step_len=self.loadtxt(self.alpha),
            func_val=self.loadtxt(self.f_try)
        )

        # New search direction needs to be searchable on disk
        if status in [0, 1]:
            m = self.load(self.m_new)
            p = self.load(self.p_new)
            self.savetxt(self.alpha, alpha)
            m_try = m + alpha * p
            self.save(self.m_try, m_try)
            self.check_model(m_try, self.m_try)

        return status

    def finalize_search(self):
        """
        Prepares algorithm machinery and scratch directory for next model update

        Removes old model/search parameters, moves current parameters to old,
        sets up new current parameters and writes statistic outputs
        """
        self.logger.info(msg.sub("FINALIZING LINE SEARCH"))

        g = self.load(self.g_new)
        p = self.load(self.p_new)
        x = self.line_search.search_history()[0]
        f = self.line_search.search_history()[1]

        # Clean scratch directory
        unix.cd(PATH.OPTIMIZE)

        # Remove the old model parameters
        if self.iter > 1:
            self.logger.info("removing previously accepted model files (old)")
            for fid in [self.m_old, self.f_old, self.g_old, self.p_old]:
                unix.rm(fid)

        self.logger.info("shifting current model (new) to previous model (old)")
        unix.mv(self.m_new, self.m_old)
        unix.mv(self.f_new, self.f_old)
        unix.mv(self.g_new, self.g_old)
        unix.mv(self.p_new, self.p_old)

        self.logger.info("setting accepted line search model as current model")
        unix.mv(self.m_try, self.m_new)
        self.savetxt(self.f_new, f.min())
        self.logger.info(f"current misfit is {self.f_new}={f.min():.3E}")

        # !!! TODO Describe what stats are being written here
        self.logger.info(f"writing optimization stats to: {CFGPATHS.STATSDIR}")
        self.write_stats(self.log_factor, value=
                         -dot(g, g) ** -0.5 * (f[1] - f[0]) / (x[1] - x[0])
                         )
        self.write_stats(self.log_gradient_norm_L1, value=np.linalg.norm(g, 1))
        self.write_stats(self.log_gradient_norm_L2, value=np.linalg.norm(g, 2))
        self.write_stats(self.log_misfit, value=f[0])
        self.write_stats(self.log_restarted, value=self.restarted)
        self.write_stats(self.log_slope, value=(f[1] - f[0]) / (x[1] - x[0]))
        self.write_stats(self.log_step_count, value=self.line_search.step_count)
        self.write_stats(self.log_step_length, value=x[f.argmin()])
        self.write_stats(self.log_theta,
                         value=180. * np.pi ** -1 * angle(p, -g))

        self.logger.info("resetting line search step count to 0")
        self.line_search.step_count = 0

    def retry_status(self):
        """
        After a failed line search, this determines if restart is worthwhile
        by checking, in effect, if the search direction was the same as gradient
        direction
        """
        g = self.load(self.g_new)
        p = self.load(self.p_new)
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
        if PAR.OPTIMIZE != "base":
            g = self.load(self.g_new)
            self.save(self.p_new, -g)
            self.line_search.clear_history()
            self.restarted = 1

    def write_stats(self, log, value=None, format="18.6E"):
        """
        Simplified write function to append values to text files in the
        STATSDIR. Used because stats line search information can be overwritten
        by subsequent iterations so we need to append values to text files
        if they should be retained.

        Log files will look something like:

        ITER  FACTOR
        ====  ======
           1     0.0

        :type log: str
        :param log: name of the file to write to. Will append .txt to it
        :type value: float
        :param value: value to write to file
        :type format: str
        :param format: string formatter for value
        """
        fid = os.path.join(PATH.WORKDIR, CFGPATHS.STATSDIR, f"{log}.txt")

        # If no value is given, assuming we are being run from setup() and
        # writing to new files. Will OVERWRITE any existing files
        if value is None:
            with open(fid, "w") as f:
                f.write(f"{'ITER':>4}  {log.upper():>18}\n")
                f.write(f"{'='*4}  {'='*18}\n")

        else:
            with open(fid, "a") as f:
                f.write(f"{self.iter:>4}  {value:{format}}\n")

    def check_model(self, m, tag):
        """
        Check to ensure that the model parameters fall within the guidelines
        of the solver. Print off min/max model parameters for the User.

        :type m: np.array
        :param m: model to check parameters of
        :type tag: str
        :param tag: tag of the model to be used for more specific error msgs
        """
        # Dynamic way to split up the model based on number of params
        pars = {}
        for i, par in enumerate(solver.parameters):
            pars[par] = np.split(m, len(solver.parameters))[i]

        # Check Poisson's ratio, which will error our SPECFEM if outside limits
        if (pars["vp"] is not None) and (pars["vs"] is not None):
            self.logger.debug(f"checking poissons ratio for: '{tag}'")
            pars["pr"] = check_poissons_ratio(vp=pars["vp"], vs=pars["vs"])
            if pars["pr"].min() < 0:
                self.logger.warning("minimum poisson's ratio is negative")

        # Tell the User min and max values of the updated model
        self.logger.info(f"model parameters ({tag} {self.eval_str}):")
        parts = "{minval:.2f} <= {key} <= {maxval:.2f}"
        for key, vals in pars.items():
            self.logger.info(parts.format(minval=vals.min(), key=key,
                                          maxval=vals.max())
                             )

    @staticmethod
    def load(filename):
        """
        Convenience function to reads vectors from disk as Numpy files,
        reads directly from PATH.OPTIMIZE. Works around Numpy's behavior of
        appending '.npy' to files that it saves.

        :type filename: str
        :param filename: filename to read from
        :rtype: np.array
        :return: vector read from disk
        """
        fid = os.path.join(PATH.OPTIMIZE, filename)
        if not os.path.exists(fid):
            fid += ".npy"
        return np.load(fid)

    @staticmethod
    def save(filename, array):
        """
        Convenience function to write vectors to disk as numpy files.
        Reads directly from PATH.OPTIMIZE

        :type filename: str
        :param filename: filename to read from
        :type array: np.array
        :param array: array to be saved
        """
        np.save(os.path.join(PATH.OPTIMIZE, filename), array)

    @staticmethod
    def loadtxt(filename):
        """
        Reads scalars from optimize directory on disk,
        accounts for savetxt() appending file extension

        :type filename: str
        :param filename: filename to read from
        :rtype: float
        :return: scalar read from disk
        """
        if not os.path.splitext(filename)[1]:
            filename += ".txt"
        return float(np.loadtxt(os.path.join(PATH.OPTIMIZE, filename)))

    @staticmethod
    def savetxt(filename, scalar):
        """
        Writes scalars to disk with a specific format, appends .txt to the
        filename to make it clear that these are text files.

        :type filename: str
        :param filename: filename to read from
        :type scalar: float
        :param scalar: value to write to disk
        """
        if not os.path.splitext(filename)[1]:
            filename += ".txt"
        np.savetxt(os.path.join(PATH.OPTIMIZE, filename), [scalar], "%11.6e")


