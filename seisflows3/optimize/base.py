#!/usr/bin/env python
"""
This is the base class for seisflows.optimize
This class provides the core utilities for the Seisflows optimization schema.
"""
import os
import sys
import logging
import numpy as np

from seisflows3.plugins import line_search, preconds
from seisflows3.tools import msg, unix
from seisflows3.tools.wrappers import loadnpy, savenpy
from seisflows3.tools.math import angle, poissons_ratio
from seisflows3.config import SeisFlowsPathsParameters, CFGPATHS


# seisflows.config objects 
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']
solver = sys.modules['seisflows_solver']


class Base:
    """
    Nonlinear optimization abstract base class

    Base class on top of which steepest descent, nonlinear conjugate, quasi-
    Newton and Newton methods can be implemented. Includes methods for both
    search direction and line search.

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
        :type restarted: int
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
        self.restarted = 0

        # Define the names of output stats logs to keep all paths in one place
        self.log_line_search = "line_search.txt"
        self.log_factor = "factor.txt"
        self.log_gradient_norm_L1 = "gradient_norm_L1.txt"
        self.log_gradient_norm_L2 = "gradient_norm_L2.txt"
        self.log_misfit = "misfit.txt"
        self.log_restarted = "restarted.txt"
        self.log_slope = "slope.txt"
        self.log_step_count = "step_count.txt"
        self.log_step_length = "step_length.txt"
        self.log_theta = "theta.txt"

        # Define the names of variables used to keep track of models etc. so
        # that we don't have multiple strings floating around defining the same
        # thing
        self.m_new = "m_new"
        self.m_old = "m_old"
        self.m_try = "m_try"
        self.f_new = "f_ew"
        self.f_old = "f_old"
        self.f_try = "f_try"
        self.g_new = "g_new"
        self.g_old = "g_old"
        self.p_new = "p_new"
        self.p_old = "p_old"
        self.alpha = "alpha"

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
                      "change in line serach behavior")

        sf.par("STEPLENINIT", required=False, default=0.05, par_type=float,
               docstr="Initial line serach step length, as a fraction "
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
        msg.check(type(self))

        if validate:
            self.required.validate()

        if PAR.OPTIMIZE == "base":
            print(msg.CompatibilityError1)
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

        # Where to write optimization statistics etc.
        path_stats = os.path.join(PATH.WORKDIR, CFGPATHS.STATSDIR)
        unix.mkdir(path_stats)

        # Prepare line search machinery
        self.line_search = getattr(line_search, PAR.LINESEARCH)(
            step_count_max=PAR.STEPCOUNTMAX,
            log_file=os.path.join(path_stats, self.log_line_search),
        )

        # Prepare preconditioner
        if PAR.PRECOND:
            self.precond = getattr(preconds, PAR.PRECOND)()
        else:
            self.precond = None

        # Ensure that line search step count starts at 0 (workflow.intialize)
        self.write_stats(self.log_step_count, 0)

        # Prepare scratch directory and save initial model
        unix.mkdir(PATH.OPTIMIZE)
        if "MODEL_INIT" in PATH:
            m_new = solver.merge(solver.load(PATH.MODEL_INIT))
            self.save(self.m_new, m_new)
            self.check_model_parameters(m_new, self.m_new)

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
        Computes search direction

        .. note::
            This function implements steepest descent, for other algorithms,
            simply overload this method
        """
        msg.whoami(type(self), prepend="computing search direction with ")

        g_new = self.load(self.g_new)
        if self.precond is not None:
            p_new = -1 * self.precond(g_new)
        else:
            p_new = -1 * g_new
        self.save(self.p_new, p_new)

    def check_model_parameters(self, m, tag):
        """
        Check to ensure that the model parameters fall within the guidelines 
        of the solver. Print off min/max model parameters for the User.

        !!! Clean this up
        
        :type m: np.array
        :param m: model to check parameters of 
        :type tag: str
        :param tag: tag of the model to be used for more specific error msgs
        """
        # Dynamic way to split up the model based on number of params
        pars = {}
        for i, par in enumerate(solver.parameters):
            pars[par] = np.split(m, len(solver.parameters))[i]

        # Check the Poisson's ratio based on Specfem3D upper/lower bounds
        if pars["vp"] is not None and pars["vs"] is not None:
            poissons = poissons_ratio(vp=pars["vp"], vs=pars["vs"])
            if (poissons.min() < -1) or (poissons.max() > 0.5):
                print(msg.PoissonsRatioError.format(tag=tag,
                                                    pmin=poissons.min(),
                                                    pmax=poissons.max())
                      )
                sys.exit(-1)
            else:
                pars["pr"] = poissons 

        # Tell the User min and max values of the updated model
        self.logger.info(f"model parameters ({tag} {self.eval_str}):")
        msg_ = "{minval:.2f} <= {key} <= {maxval:.2f}"
        for key, vals in pars.items():
            self.logger.info(msg_.format(minval=vals.min(), key=key,
                                         maxval=vals.max()))

    def initialize_search(self):
        """
        Determines first step length in line search
        """
        # Load in and calucate the necessary variables
        m = self.load(self.m_new)
        g = self.load(self.g_new)
        p = self.load(self.p_new)
        f = self.loadtxt(self.f_new)
        norm_m = max(abs(m))
        norm_p = max(abs(p))
        gtg = self.dot(g, g)
        gtp = self.dot(g, p)

        # Restart line search if necessary
        if self.restarted:
            self.line_search.clear_history()

        # Optional step length safeguard
        if PAR.STEPLENMAX:
            self.line_search.step_len_max = PAR.STEPLENMAX * norm_m / norm_p

        # Determine initial step length
        alpha, _ = self.line_search.initialize(iter=self.iter, step_len=0.,
                                               func_val=f, gtg=gtg, gtp=gtp
                                               )

        # Optional initial step length override
        if PAR.STEPLENINIT and len(self.line_search.step_lens) <= 1:
            alpha = PAR.STEPLENINIT * norm_m / norm_p
            self.logger.debug(f"step length override due to "
                              f"PAR.STEPLENINIT={PAR.STEPLENINIT}")

        # The new model is the old model, scaled by the step direction and
        # gradient threshold to remove any outlier values
        m_try = m + alpha * p

        # Write model corresponding to chosen step length
        self.save(self.m_try, m_try)
        self.savetxt(self.alpha, alpha)

        # Check the new model and update the User on a few parameters
        self.check_model_parameters(m_try, self.m_try)

    def update_search(self):
        """
        Updates line search status and step length

        Status codes:
            status > 0  : finished
            status == 0 : not finished
            status < 0  : failed
        """
        alpha, status = self.line_search.update(
            iter=self.iter, step_len=self.loadtxt(self.alpha),
            func_val=self.loadtxt(self.f_try)
        )

        # write model corresponding to chosen step length
        if status >= 0:
            m = self.load(self.m_new)
            p = self.load(self.p_new)
            self.savetxt(self.alpha, alpha)
            m_try = m + alpha * p
            self.save(self.m_try, m_try)
            self.check_model_parameters(m_try, self.m_try)

        return status

    def finalize_search(self):
        """
        Prepares algorithm machinery and scratch directory for next model update

        Removes old model/search parameters, moves current parameters to old,
        sets up new current parameters and writes statistic outputs
        """
        # m = self.load('m_new')  # unusued variable
        g = self.load(self.g_new)
        p = self.load(self.p_new)
        x = self.line_search.search_history()[0]
        f = self.line_search.search_history()[1]

        # Clean scratch directory
        unix.cd(PATH.OPTIMIZE)
        # Remove the old model parameters
        if self.iter > 1:
            for fid in [self.m_old, self.f_old, self.g_old, self.p_old]:
                unix.rm(fid)

        # Rename current model parameters to "_old" for new search
        unix.mv(self.m_new, self.m_old)
        unix.mv(self.f_new, self.f_old)
        unix.mv(self.g_new, self.g_old)
        unix.mv(self.p_new, self.p_old)

        # Setup the current model parameters
        unix.mv(self.m_try, self.m_new)
        self.savetxt(self.f_new, f.min())

        # Output the latest statistics to text files
        # !!! Describe what stats are being written here
        self.write_stats(filename=self.log_factor, value=
                         -self.dot(g, g) ** -0.5 * (f[1] - f[0]) / (x[1] - x[0])
                         )
        self.write_stats(filename=self.log_gradient_norm_L1,
                         value=np.linalg.norm(g, 1))
        self.write_stats(filename=self.log_gradient_norm_L2,
                         value=np.linalg.norm(g, 2))
        self.write_stats(filename=self.log_misfit, value=f[0])
        self.write_stats(filename=self.log_restarted, value=self.restarted)
        self.write_stats(filename=self.log_slope,
                         value=(f[1] - f[0]) / (x[1] - x[0]))
        self.write_stats(filename=self.log_step_count,
                         value=self.line_search.step_count)
        self.write_stats(filename=self.log_step_length,
                         value=x[f.argmin()])
        self.write_stats(filename=self.log_theta,
                         value=180. * np.pi ** -1 * angle(p, -g))

        # Reset line search step count to 0 for next iteration
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
        Restarts nonlinear optimization algorithm.
        Keeps current position in model space, but discards history of
        nonlinear optimization algorithm in an attempt to recover from
        numerical stagnation.
        """
        g = self.load(self.g_new)
        self.save(self.p_new, -g)
        self.line_search.clear_history()
        self.restarted = 1

    def write_stats(self, filename, value, format="e"):
        """
        Simplified write function to append values to text files in the
        STATSDIR. Used because stats line search information can be overwritten
        by subsequent iterations so we need to append values to text files
        if they should be retained.

        :type filename: str
        :param filename: name of the file to write to
        :type value: float
        :param value: value to write to file
        :type format: str
        :param format: string formatter for value
        """
        fid = os.path.join(CFGPATHS.STATSDIR, {filename})
        if not os.path.exists(fid):
            self.logger.debug(f"creating stats file: {fid}")
        with open(fid, "a") as f:
            f.write(f"{value:{format}}\n")

    @staticmethod
    def dot(x, y):
        """
        Utility function to computes inner product between vectors

        :type x: np.array
        :param x: vector 1
        :type y: np.array
        :param y: vector 2
        """
        return np.dot(np.squeeze(x), np.squeeze(y))

    @staticmethod
    def load(filename):
        """
        Reads vectors from disk

        :type filename: str
        :param filename: filename to read from
        :return:
        """
        fid_out = os.path.join(PATH.OPTIMIZE, filename)
        # logger.debug(f"loading model vector from: {fid_out}")
        return loadnpy(os.path.join(PATH.OPTIMIZE, filename))

    @staticmethod
    def save(filename, array):
        """
        Writes vectors to disk

        :type filename: str
        :param filename: filename to read from
        :type array: np.array
        :param array: array to be saved
        :return:
        """
        fid_out = os.path.join(PATH.OPTIMIZE, filename)
        # logger.debug(f"saving model vector to: {fid_out}")
        savenpy(fid_out, array)

    @staticmethod
    def loadtxt(filename):
        """
        Reads scalars from disk

        :type filename: str
        :param filename: filename to read from
        :return:
        """
        return float(np.loadtxt(os.path.join(PATH.OPTIMIZE, filename)))

    @staticmethod
    def savetxt(filename, scalar):
        """
        Writes scalars to disk

        :type filename: str
        :param filename: filename to read from
        :type scalar: float
        :param scalar: value to write to disk
        :return:
        """
        np.savetxt(os.path.join(PATH.OPTIMIZE, filename), [scalar], "%11.6e")


