#!/usr/bin/env python3
"""
Gradient descent nonlinear optimization algorithm. Acts as the Base class for
optimization.

The Optimization library contains classes and methods used to solve nonlinear
optimization problems, i.e., misfit minimization. Various subclasses implement
different optimization algorithms.

.. note::
    To reduce memory overhead, vectors are read from disk rather than passed
    from calling routines. For example, at the beginning of
    `compute_direction` the current gradient is read from  'g_new' and the
    resulting search direction is written to 'p_new'. As the inversion
    progresses, other information is stored as well.

.. note::
    The default numerical parameters for each algorithm should work well for a
    range of applications without manual tuning. If the nonlinear
    optimization procedure stagnates, it may be due to issues involving
    data quality or the choice of data misfit, data processing, or
    regularization parameters.  Problems in any of these areas usually
    manifest themselves through stagnation of the nonlinear optimization
    algorithm.

TODO fix line search, currently broken. Store line search variables in .npz
    arrays rather than as internal attributes, that way it's easier to restart
    from a broken workflow. Need to figure out how to reset and restart
    line searches efficiently
"""
import os
import sys
import numpy as np

from seisflows import logger
from seisflows.tools import msg, unix
from seisflows.tools.core import Dict
from seisflows.tools.math import angle, dot
from seisflows.tools.specfem import Model
from seisflows.plugins import line_search as line_search_dir


class Gradient:
    """
    [optimization.gradient] gradient/steepest descent optimization algorithm.

    :type line_search_method: str
    :param line_search_method: chosen line_search algorithm. Currently available
        are 'bracket' and 'backtrack'. See seisflows.plugins.line_search
        for all available options
    :type preconditioner: str
    :param preconditioner: algorithm for preconditioning gradients
    :type step_count_max: int
    :param step_count_max: maximum number of trial steps to perform during
        the line search before a change in line search behavior is
        considered
    :type step_len_init: float
    :param step_len_init: initial line search step length as a fraction of
        current model parameters.
    :type step_len_max: float
    :param step_len_max: maximum allowable step length during the line
        search. Set as a fraction of the current model parameters
    """
    def __init__(self, start=1, line_search_method="bracket", preconditioner=None,
                 step_count_max=10, step_len_init=0.05, step_len_max=0.5,
                 workdir=os.getcwd(), path_optimize=None, path_output=None,
                 path_preconditioner=None,
                 **kwargs):
        """Gradient-descent input parameters"""
        super().__init__()

        self.iteration = start  # to match PAR.BEGIN
        self.preconditioner = preconditioner

        self.step_count_max = step_count_max
        self.step_len_init = step_len_init
        self.step_len_max = step_len_max

        # Internal check to see if the chosen line search algorithm exists
        if not hasattr(line_search_dir, line_search_method):
            logger.warning(f"{line_search_method} is not a valid line search "
                           f"algorithm, defaulting to 'bracket'")
            line_search_method = "bracket"

        # .title() ensures we grab the class and not the module
        self.line_search = getattr(line_search_dir, line_search_method.title())(
            step_count_max=step_count_max, step_len_max=step_len_max
        )

        # Set required path structure
        self.path = Dict(
            scratch=path_optimize or
                    os.path.join(os.getcwd(), "scratch", "optimize"),
            output=path_output or os.path.join(workdir, "output"),
            preconditioner=path_preconditioner,

        )

        # Internally used parameters for checking validity
        self._acceptable_vectors = ["m_new", "m_old", "m_try",
                                    "g_new", "g_old", "g_try",
                                    "p_new", "p_old", "alpha",
                                    "f_new", "f_old", "f_try"]
        self._acceptable_preconditioners = ["diagonal"]

        self.restarted = False

    @property
    def step_count(self):
        """Convenience property to access `step_count` from line search"""
        return self.line_search.step_count

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

        assert 0. < self.step_len_init, f"optimize.step_len_init must be >= 0."
        assert 0. < self.step_len_max, f"optimize.step_len_max must be >= 0."
        assert self.step_len_init < self.step_len_max, \
            f"optimize.step_len_init must be < optimize.step_len_max"

    def setup(self):
        """
        Sets up nonlinear optimization machinery
        """
        unix.mkdir(self.path.scratch)

    def load(self, name):
        """
        Convenience function to access the full paths of model and gradient
        vectors that are saved to disk

        .. note:: the available options that can be loaded
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
        model_txt = model_npz.replace(".npz", ".txt")
        if os.path.exists(model_npz):
            model = Model(path=os.path.join(self.path.scratch,
                                            f"{name}.npz"), load=True)
        elif os.path.exists(model_txt):
            model = float(np.loadtxt(model_txt))
        else:
            raise FileNotFoundError(f"no optimization file found for {name}")

        return model

    def save(self, name, m):
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
        elif isinstance(m, (float, int)):
            path = os.path.join(self.path.scratch, f"{name}.txt")
            np.savetxt(path, [m])
        else:
            raise TypeError(f"optimize.save unrecognized type error {type(m)}")

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
        """
        g_new = self.load("g_new")
        p_new = g_new.update(vector=-1 * self._precondition(g_new.vector))

        return p_new

    def initialize_search(self):
        """
        Initialize the plugin line search machinery. Should only be run at
        the beginning of line search, by the main workflow module.
        """
        # Vectors required to initialize a line search
        m = self.load("m_new")
        g = self.load("g_new")
        p = self.load("p_new")
        f = self.load("f_new")

        norm_m = max(abs(m.vector))
        norm_p = max(abs(p.vector))
        gtg = dot(g.vector, g.vector)
        gtp = dot(g.vector, p.vector)

        # Restart plugin line search if the optimization library restarts
        if self.restarted:
            self.line_search.clear_history()

        # Optional safeguard to prevent step length from getting too large
        if self.step_len_max:
            new_step_len_max = self.step_len_max * norm_m / norm_p
            self.line_search.step_len_max = new_step_len_max
            logger.debug(f"max step length safeguard = {new_step_len_max:.2E}")

        # Initialize the line search by setting vector variables
        self.line_search.step_count = 0
        self.line_search.step_lens.append(0.)
        self.line_search.func_vals.append(f)
        self.line_search.gtg.append(gtg)
        self.line_search.gtp.append(gtp)

        alpha, _ = self.line_search.calculate_step()

        # Alpha defines the trial step length. Optional step length override
        if self.step_len_init and len(self.line_search.step_lens) <= 1:
            alpha = self.step_len_init * norm_m / norm_p
            logger.debug(f"overwrite initial step length: {alpha:.2E}")

        # The new model is the old model, scaled by the step direction and
        # gradient threshold to remove any outlier values
        m_try = m.update(vector=m.vector + alpha * p.vector)

        return m_try, alpha

    def update_search(self):
        """
        Updates line search status and step length and checks if the line search
        has been completed.

        Available status codes from line_search.calculate_step():
            status == 1  : finished
            status == 0 : not finished
            status == -1  : failed
        """
        # Update the line search variables for a new step calculation
        self.line_search.step_lens.append(self.load("alpha"))
        self.line_search.func_vals.append(self.load("f_try"))

        alpha, status = self.line_search.calculate_step()

        # New search direction needs to be searchable on disk
        if status in [0, 1]:
            m = self.load("m_new")
            p = self.load("p_new")
            self.save("alpha", alpha)

            m_try = m.update(vector=m.vector + alpha * p.vector)
            self.save("m_try", m_try)

        return status

    def finalize_search(self):
        """
        Prepares algorithm machinery and scratch directory for next model update

        Removes old model/search parameters, moves current parameters to old,
        sets up new current parameters and writes statistic outputs
        """
        unix.cd(self.path.scratch)

        logger.info(msg.sub("FINALIZING LINE SEARCH"))

        # Remove the old model parameters
        if self.iteration > 1:
            logger.info("removing previously accepted model files (old)")
            for fid in ["m_old", "f_old", "g_old", "p_old"]:
                unix.rm(fid)

        # Needs to be run before shifting model in next step
        self._write_stats()

        logger.info("shifting current model (new) to previous model (old)")
        unix.mv("m_new.npz", "m_old.npz")
        unix.mv("f_new.npz", "f_old.npz")
        unix.mv("g_new.npz", "g_old.npz")
        unix.mv("p_new.npz", "p_old.npz")

        logger.info("setting accepted line search model as current model")
        unix.mv("m_try.npz", "m_new.npz")

        # Choose minimum misfit value as final misfit/model
        f = self.line_search.get_search_history()[1]
        self.save("f_new", f.min())
        logger.info(f"current misfit is {f.min():.3E}")

        logger.info("resetting line search step count to 0")
        self.line_search.step_count = 0

    def retry_status(self):
        """
        After a failed line search, this determines if restart is worthwhile
        by checking, in effect, if the search direction was the same as gradient
        direction
        """
        g = self.load("g_new")
        p = self.load("p_new")
        theta = angle(p.vector, -1 * g.vector)

        logger.debug(f"theta: {theta:6.3f}")

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

        .. note::
            steepest descent optimization algorithm does not have any restart
            capabilities. This function is instantiated here to be overwritten
            by child classes
        """
        pass

    def _write_stats(self):
        """
        Simplified write function to append values to text files.
        Used because stats line search information can be overwritten
        by subsequent iterations so we need to append values to text files
        if they should be retained.
        """
        logger.info(f"writing optimization stats")
        fid = os.path.join(self.path.output,  f"optim_stats.txt")

        # First time, write header information
        if not os.path.exists(fid):
            with open(fid, "w") as f:
                for header in ["ITER", "FACTOR", "GRAD_NORM_L1", "GRAD_NORM_L2",
                               "MISFIT", "RESTART", "SLOPE", "STEP", "LENGTH",
                               "THETA"]:
                    f.write(f"{header.upper()},")
                f.write("\n")

        g = self.load("g_new")
        p = self.load("p_new")
        x, f, *_ = self.line_search.get_search_history()

        import pdb;pdb.set_trace()

        # Calculated stats factors
        factor = -1 * dot(g.vector, g.vector)
        factor = factor ** -0.5 * (f[1] - f[0]) / (x[1] - x[0])

        grad_norm_L1 = np.linalg.norm(g.vector, 1)
        grad_norm_L2 = np.linalg.norm(g.vector, 2)

        misfit = f[0]
        restarted = self.restarted
        slope = (f[1] - f[0]) / (x[1] - x[0])
        step_count = self.line_search.step_count
        step_length = x[f.argmin()]
        theta = 180. * np.pi ** -1 * angle(p, -1 * g.vector)

        with open(fid, "a") as f:
            f.write(f"{self.iteration:0>2},"
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
