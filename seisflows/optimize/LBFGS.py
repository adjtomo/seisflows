#!/usr/bin/env python3
"""
This is the custom class for an LBFGS optimization schema.
It supercedes the `seisflows.optimize.base` class
"""
import os
import numpy as np

from seisflows import logger
from seisflows.optimize.gradient import Gradient
from seisflows.tools import unix
from seisflows.tools.msg import DEG
from seisflows.tools.math import angle
from seisflows.plugins import line_search as line_search_dir


class LBFGS(Gradient):
    """
    The Limited memory BFGS algorithm
    Calls upon seisflows.plugin.optimize.LBFGS to accomplish LBFGS algorithm

    Includes optional safeguards: periodic restarting and descent conditions.

    To conserve memory, most vectors are read from disk rather than passed
    from a calling routine.

    L-BFGS Variables:
        s: memory of model differences
        y: memory of gradient differences

    Optimization Variables:
        m: model
        f: objective function value
        g: gradient direction
        p: search direction

    Line Search Variables:
        x: list of step lenths from current line search
        f: correpsonding list of function values
        m: number of step lengths in current line search
        n: number of model updates in optimization problem
        gtg: dot product of gradient with itself
        gtp: dot product of gradient and search direction

    Status codes
        status > 0  : finished
        status == 0 : not finished
        status < 0  : failed
    """
    def __init__(self, lbfgs_mem=3, lbfgs_max=np.inf, lbfgs_thresh=0.,
                 **kwargs):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type lbfgs_mem: int
        :param lbfgs_mem: L-BFGS memory. Max number of previous gradients to
            retain in local memory for approximating the objective function.
        :type lbfgs_max: L-BFGS periodic restart interval. Must be
            1 <= lbfgs_max <= infinity.
        :type lbfgs_thresh: L-BFGS angle restart threshold. If the angle between
            the current and previous search direction exceeds this value,
            optimization algorithm will be restarted.


        :type LBFGS_iter: int
        :param LBFGS_iter: an internally used iteration that differs from
            optimization iter. Keeps track of internal LBFGS memory of previous
            gradients. If LBFGS is restarted, the LBFGS_iter iteration is reset,
            but the optization iteration.
        :type memory_used: int
        :param memory_used: bookkeeping to see how many previous
            gradients have been stored to internal memory. Should not exceed
            PAR.LBFGSMEM
        :type LBFGS_dir: str
        :param LBFGS_dir: location to store LBFGS internal memory
        :type y_file: str
        :param y_file: path to store memory of the gradient differences
            i.e., `g_new - g_old`
        :type s_file: str
        :param s_file: path to store memory of the model differences
            i.e., `m_new - m_old`
        """
        super().__init__(**kwargs)

        # Overwrite user-chosen line search. L-BFGS requires 'Backtrack'ing LS
        if self._line_search.title != "Backtrack":
            logger.warning(f"L-BFGS optimization requires 'backtrack'ing line "
                           f"search. Overwritng {self._line_search}")
            self._line_search = "Backtrack"
            self.line_search = getattr(line_search_dir, self._line_search)(
                step_count_max=self.step_count_max,
                step_len_max=self.step_len_max
            )

        self.LBFGS_mem = lbfgs_mem
        self.LBFGS_max = lbfgs_max
        self.LBFGS_thresh = lbfgs_thresh

        # Internally used memory and path parameters
        self._LBFGS_iter = 0
        self._memory_used = 0
        self._LBFGS_dir = os.path.join(self.path, "LBFGS")
        self._y_file = os.path.join(self.path, "LBFGS", "Y")
        self._s_file = os.path.join(self.path, "LBFGS", "S")

    def setup(self):
        """
        Set up the LBFGS optimization schema
        """
        super().setup()
        unix.mkdir(self._LBFGS_dir)

    def compute_direction(self):
        """
        Call on the L-BFGS optimization machinery to compute a search
        direction using internally stored memory of previous gradients.
        The potential outcomes when computing direction with L-BFGS

        TODO do we need to precondition L-BFGS?

        1. First iteration of L-BFGS optimization, search direction is defined
            as the inverse gradient
        2. L-BFGS internal iteration ticks over the maximum allowable number of
            iterations, force a restart condition, search direction is the
            inverse gradient
        3. New search direction vector is too far from previous direction,
            force a restart, search direction is inverse gradient
        4. New search direction is acceptably angled from previous,
            becomes the new search direction
        """
        unix.cd(self.path)

        logger.info(f"computing search direction with L-BFGS")
        self._LBFGS_iter += 1

        # Load the current gradient direction, which is the L-BFGS search
        # direction if this is the first iteration
        g = self.load("g_new")
        if self._LBFGS_iter == 1:
            logger.info("first L-BFGS iteration, setting search direction "
                             "as inverse gradient")
            p_new = -1 * g.vector
            restarted = 0

        # Restart condition or first iteration lead to setting search direction
        # as the inverse gradient (i.e., default to steepest descent)
        elif self._LBFGS_iter > self.LBFGS_max:
            logger.info("restarting L-BFGS due to periodic restart condition. "
                        "setting search direction as inverse gradient")
            self.restart()
            p_new = -1 * g.vector
            restarted = 1
        # Normal LBFGS direction computation
        else:
            # Update the search direction, apply the inverse Hessian such that
            # 'q' becomes the new search direction 'g'
            logger.info("applying inverse Hessian to gradient")
            s, y = self._update()
            q = self._apply(g.vector, s, y)

            # Determine if the new search direction is appropriate by checking
            # its angle to the previous search direction
            if self._check_status(g, q):
                logger.info("new L-BFGS search direction found")
                p_new = -q
                restarted = 0
            else:
                logger.info("new search direction not appropriate, defaulting "
                            "to inverse gradient")
                self.restart()
                p_new = -g
                restarted = 1

        # Save values to disk and memory
        self.save("p_new", p_new)
        self.restarted = restarted

    def restart(self):
        """
        Overwrite the optimization restart to restart the L-BFGS schema
        """
        super().restart()

        logger.info("restarting L-BFGS optimization algorithm by clearing "
                         "internal memory")
        self._LBFGS_iter = 1
        self._memory_used = 0

        unix.cd(self.path.OPTIMIZE)
        s = np.memmap(filename=self._s_file, mode="r+")
        y = np.memmap(filename=self._y_file, mode="r+")
        s[:] = 0.
        y[:] = 0.

    def _update(self):
        """
        Updates L-BFGS algorithm history

        .. note::
            Because models are large, and multiple iterations of models need to
            be stored in memory, previous models are stored as `memmaps`,
            which allow for access of small segments of large files on disk,
            without reading the entire file. Memmaps are array like objects.

        .. note::
            Notation for s and y taken from Liu & Nocedal 1989
            iterate notation: sk = x_k+1 - x_k and yk = g_k+1 - gk

        :rtype s: np.memmap
        :return s: memory of the model differences `m_new - m_old`
        :rtype y: np.memmap
        :return y: memory of the gradient differences `g_new - g_old`
        """
        unix.cd(self.path)

        # Determine the iterates for model m and gradient g
        s_k = self.load("m_new").vector - self.load("m_old").vector
        y_k = self.load("g_new").vector - self.load("g_old").vector

        # Determine the shape of the memory map (length of model, length of mem)
        m = len(s_k)
        n = self.LBFGS_mem

        # Initial iteration, need to create the memory map
        if self._memory_used == 0:
            s = np.memmap(filename=self._s_file, mode="w+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename=self._y_file, mode="w+", dtype="float32",
                          shape=(m, n))
            # Store the model and gradient differences in memmaps
            s[:, 0] = s_k
            y[:, 0] = y_k
            self._memory_used = 1
        # Subsequent iterations will append to memory maps
        else:
            s = np.memmap(filename=self._s_file, mode="r+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename=self._y_file, mode="r+", dtype="float32",
                          shape=(m, n))
            # Shift all stored memory by one index to make room for latest mem
            s[:, 1:] = s[:, :-1]
            y[:, 1:] = y[:, :-1]
            # Store the latest model and gradient in first index
            s[:, 0] = s_k
            y[:, 0] = y_k

            # Keep track of the memory used
            if self._memory_used < self.LBFGS_mem:
                self._memory_used += 1

        return s, y

    def _apply(self, q, s=None, y=None):
        """
        Applies L-BFGS inverse Hessian to given vector

        :type q: np.array
        :param q: gradient direction to apply L-BFGS to
        :type s: np.memmap
        :param s: memory of model differences
        :param s: memory of model differences
        :type y: np.memmap
        :param y: memory of gradient direction differences
        :rtype r: np.array
        :return r: new search direction from application of L-BFGS
        """
        unix.cd(self.path)

        # If no memmaps are given as arguments, instantiate them
        if s is None or y is None:
            m = len(q)
            n = self.LBFGS_mem
            s = np.memmap(filename=self._s_file, mode="w+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename=self._y_file, mode="w+", dtype="float32",
                          shape=(m, n))

        # First matrix product
        # Recursion step 2 from appendix A of Modrak & Tromp 2016
        kk = self._memory_used
        rh = np.zeros(kk)
        al = np.zeros(kk)
        for ii in range(kk):
            rh[ii] = 1 / np.dot(y[:, ii], s[:, ii])
            al[ii] = rh[ii] * np.dot(s[:, ii], q)
            q = q - al[ii] * y[:, ii]

        # Apply a preconditioner if available
        if self.preconditioner:
            r = self._precondition(q)
        else:
            r = q

        # Use scaling M3 proposed by Liu and Nocedal 1989
        sty = np.dot(y[:, 0], s[:, 0])
        yty = np.dot(y[:, 0], y[:, 0])
        r *= sty/yty

        # Second matrix product
        # Recursion step 4 from appendix A of Modrak & Tromp 2016
        for ii in range(kk - 1, -1, -1):
            be = rh[ii] * np.dot(y[:, ii], r)
            r = r + s[:, ii] * (al[ii] - be)

        return r

    def _check_status(self, g, r):
        """
        Check the status of the apply() function, determine if restart necessary
        Return of False means restart, return of True means good to go.

        :type g: np.array
        :param g: current gradient direction
        :type r: np.array
        :param r: new gradient direction
        :rtype: bool
        :return: okay status based on status check (False==bad, True==good)
        """
        theta = 180. * np.pi ** -1 * angle(g, r)
        logger.info(f"new search direction: {theta:.2f}{DEG} from current")

        if not 0. < theta < 90.:
            logger.info("restarting L-BFGS, theta not a descent direction")
            okay = False
        elif theta > 90. - self.LBFGS_thresh:
            logger.info("restarting L-BFGS due to practical safeguard")
            okay = False
        else:
            okay = True
        return okay
