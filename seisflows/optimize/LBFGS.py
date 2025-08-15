#!/usr/bin/env python3
"""
L-BFGS (Limited memory Broyden–Fletcher–Goldfarb–Shanno) algorithm for solving
nonlinear optimization problems.

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

TODO store LBFGS_iter during checkpointing
"""
import os
import numpy as np

from seisflows import logger
from seisflows.optimize.gradient import Gradient
from seisflows.tools import unix
from seisflows.tools.msg import DEG
from seisflows.tools.math import angle
from seisflows.tools.specfem_model import Model
from seisflows.plugins import line_search as line_search_dir


class LBFGS(Gradient):
    """
    L-BFGS Optimization
    -------------------
    Limited memory BFGS nonlinear optimization algorithm

    Parameters
    ----------
    :type lbfgs_mem: int
    :param lbfgs_mem: L-BFGS memory. Max number of previous gradients to
        retain in local memory for approximating the objective function.
    :type lbfgs_max: L-BFGS periodic restart interval. Must be
        1 <= lbfgs_max <= infinity.
    :type lbfgs_thresh: L-BFGS angle restart threshold. If the angle between
        the current and previous search direction exceeds this value,
        optimization algorithm will be restarted.

    Paths
    -----
    ***
    """
    __doc__ = Gradient.__doc__ + __doc__

    def __init__(self, lbfgs_mem=3, lbfgs_max=np.inf, lbfgs_thresh=0.,
                 **kwargs):
        """Instantiate L-BFGS specific parameters"""
        super().__init__(**kwargs)

        # Overwrite user-chosen line search. L-BFGS requires 'Backtrack'ing LS
        if self.line_search_method.title() != "Backtrack":
            self.line_search_method = "Backtrack"
            self._line_search = getattr(
                line_search_dir, self.line_search_method)(
                    step_count_max=self.step_count_max,
            )

        self.LBFGS_mem = lbfgs_mem
        self.LBFGS_max = lbfgs_max
        self.LBFGS_thresh = lbfgs_thresh

        # Set new L-BFGS dependent paths for storing previous gradients
        self.path["_LBFGS"] = os.path.join(self.path.scratch, "LBFGS")
        # self.path["_y_file"] = os.path.join(self.path["_LBFGS"], "Y.dat")
        # self.path["_s_file"] = os.path.join(self.path["_LBFGS"], "S.dat")

        # These are used to store LBFGS memory which include model and gradient
        # differences
        for i in self.LBFGS_mem:
            y = f"_y{i}"  # `y` = model difference (m_i+1 - m_i)
            self.path[y] = os.path.join(self.path._LBFGS, y)

            s = f"_s{i}"  # `s` = gradient difference (g_i+1 - g_i)
            self.path[s] = os.path.join(self.path._LBFGS, s)

        # Internally used memory parameters for the L-BFGS optimization algo.
        self._LBFGS_iter = 0
        self._memory_used = 0

    def setup(self):
        """
        Set up the LBFGS optimization schema
        """
        super().setup()
        unix.mkdir(self.path._LBFGS)

    def checkpoint(self):
        """
        Overwrite default checkpointing to store internal L-BFGS Attributes
        """
        super().checkpoint()
        checkpoint_dict = dict(np.load(self.path._checkpoint))
        checkpoint_dict["LBFGS_iter"] = self._LBFGS_iter
        checkpoint_dict["memory_used"] = self._memory_used

        np.savez(file=self.path._checkpoint, **checkpoint_dict)  # NOQA

    def load_checkpoint(self):
        """
        Counterpart to `optimize.checkpoint`. Loads a checkpointed optimization
        module from disk and sets internal tracking attributes. Adds additional
        functionality to restore internal L-BFGS attributes
        """
        super().load_checkpoint()

        # NumPy appends '.npz' when saving. Make sure we honor that.
        if not self.path._checkpoint.endswith(".npz"):
            fid = f"{self.path._checkpoint}.npz"
        else:
            fid = self.path._checkpoint

        if os.path.exists(fid):
            dict_in = np.load(file=fid)
            self._LBFGS_iter = int(dict_in["LBFGS_iter"])
            self._memory_used = int(dict_in["memory_used"])

    def compute_direction(self):
        """
        Call on the L-BFGS optimization machinery to compute a search
        direction using internally stored memory of previous gradients.

        The potential outcomes when computing direction with L-BFGS:
        1. First iteration of L-BFGS optimization, search direction is defined
            as the inverse gradient
        2. L-BFGS internal iteration ticks over the maximum allowable number of
            iterations, force a restart condition, search direction is the
            inverse gradient
        3. New search direction vector is too far from previous direction,
            force a restart, search direction is inverse gradient
        4. New search direction is acceptably angled from previous,
            becomes the new search direction

        :rtype: seisflows.tools.specfem.Model
        :return: search direction as a Model instance
        """
        self._LBFGS_iter += 1
        logger.info(f"LBFGS iteration incremented -> {self._LBFGS_iter}")

        # Load the current gradient direction, which is the L-BFGS search
        # direction if this is the first iteration
        g = Model(path=self.path._g_new)

        if self._LBFGS_iter == 1:
            logger.info("first L-BFGS iteration, default to gradient descent "
                        "(P = -G)")
            # Search direction `p_new` is negative gradient 
            g.apply(actions=["*"], values=[-1], export_to=self.path._p_new)
            restarted = False
        # Restart condition or first iteration lead to setting search direction
        # as the inverse gradient (i.e., default to steepest descent)
        elif self._LBFGS_iter > self.LBFGS_max:
            logger.info("restarting L-BFGS due to periodic restart condition. "
                        "setting search direction as inverse gradient (P = -G)")
            self.restart()
            restarted = True

            # New search direction `p_new` is negative gradient
            g.apply(actions=["*"], values=[-1], export_to=self.path._p_new)
        # Normal LBFGS direction computation
        else:
            # Update the search direction, apply the inverse Hessian such that
            # 'q' becomes the new search direction 'g'
            logger.info("applying inverse Hessian to gradient")
            self.update_search_history()
            q = g.copy()
            q. update(vector=self._apply_inverse_hessian(g.vector, s, y))

            # Determine if the new search direction is appropriate by checking
            # its angle to the previous search direction
            if self._check_status(g.vector, q.vector):
                logger.info("new L-BFGS search direction found")
                p_new.update(vector=-1 * q.vector)
                restarted = False
            else:
                logger.info("new search direction not appropriate, defaulting "
                            "to gradient descent")
                self.restart()
                p_new.update(vector=-1 * g.vector)
                restarted = True

        # Assign restart condition to internal memory
        self._restarted = restarted

        return p_new

    def restart(self):
        """
        Restart the L-BFGS optimization algorithm by clearing out stored
        gradient memory. The workflow will need to call `compute_direction` and
        `initialize_search` afterwards to properly re-instantiate the line
        search machinery. 
        """
        logger.info("restarting L-BFGS optimization algorithm")

        # Clear internal memory back to default values
        self._line_search.clear_search_history()
        self._restarted = True
        self._LBFGS_iter = 0  # will be incremented to 1 by `compute_direction`
        self._memory_used = 0

        # Clear out previous gradient information. Check file existence first
        # because the first iteration will not have generated S and Y yet, but
        # we may want to run restart for manual line search restart
        if os.path.exists(self.path._s_file):
            s = np.memmap(filename=self.path._s_file, mode="r+")
            s[:] = 0.
        if os.path.exists(self.path._y_file):
            y = np.memmap(filename=self.path._y_file, mode="r+")
            y[:] = 0.

    def update_search_history(self):
        """
        Updates L-BFGS algorithm history by removing old memory and 

        .. note::
            This algorithm is notated in Modrak & Tromp 2016 Appendix A
        
        .. note::
            Notation for s and y taken from Liu & Nocedal 1989
            iterate notation: sk = x_k+1 - x_k and yk = g_k+1 - gk
        """
        # If we have memory assigned, we need to make space for the latest entry
        # Essentially were shifting every index down by one (i -> i+1) and 
        # removing the last index to make way for new 0th index memory. In 
        # practice we're just renaming directories
        if self._memory_used > 0:
            if self._memory_used == self.LBFGS_mem:
                # Remove the oldest memory if we exceed the memory threshold
                unix.rm(self.path._LBFGS[f"_y{self._memory_used - 1}"])
                unix.rm(self.path._LBFGS[f"_s{self._memory_used - 1}"])

            # Iterate backwards, so for mem==3 we go: 2, 1 (ignore 0)
            for i in np.arange(self._memory_used - 1, 0, -1):
                # Reassigning the next index to the current. E.g., if mem==3 we
                # start at i==2 which we don't need anymore, and move i==1 into 
                # the 2 spot.
                unix.rename(old=self.path.LBFGS[f"_y{i-1}"], 
                            new=self.path.LBFGS[f"_y{i}"])
                unix.rename(old=self.path.LBFGS[f"_s{i-1}"],
                            new=self.path.LBFGS[f"_s{i}"])
        
        # Assign the current model and gradient differences to 0th index
        # s_k = m_new - m_old
        Model(self.path._m_new).apply(
            actions=["-"], other=Model(self.path._m_old),
            export_to=self.path.LBFGS._s0
            )
        # y_k = g_new - g_old
        Model(self.path._g_new).apply(
            actions=["-"], other=Model(self.path._g_old),
            export_to=self.path.LBFGS._y0
            )

        # Increment memory up to LBFGS memory threshold
        if self._memory_used < self.LBFGS_mem:
            self._memory_used += 1

    def apply_inverse_hessian(self):
        """
        Applies L-BFGS inverse Hessian to given vector following Appendix A
        Recursion from Modrak and Tromp (2016)

        In words the Recursion step is as follows:
        (1) Set q = gk
            i = k− 1
            j= min (k, l) 
        (2) Perform j times: 
            lambda_i = s_i^T * q / y_i^T * s_i, 
            q -= lambda_i * y_i
            i -= 1
        (3) Set gamma = s^T_(k-1) * y_(k-1) / y^T_(k-1) * y_k-1
            r = gamma * D * q
            i = k - j
        (4) Perform j times
            mu = y_i^T * r / y_i^T
            r += s_i(lambda_i - mu)
            i += 1
        (5) End with result: p_k = -r
        """
        # Step 1: Set q=g_k, i=k-1, j=min(k,l)
        q = Model(path=self.path._g_new)
        k = self._memory_used
        j = min(k, self.LBFGS_mem)
        
        # Step 2: Perform j times
        for i in range(j):
            y_i = Model(self.path.LBFGS[f"_y_{i}"])
            s_i = Model(self.path.LBFGS[f"_s_{i}"])

            # lambda_i = s_i^T * q / y_i^T * s_i
            lambda_i = s_i.dot(q) / y_i.dot(s_i)

            # q -= lambda_i * y_i
            q.apply(actions=["-"], values=[al])



    # def _update_search_history(self):
    #     """
    #     Updates L-BFGS algorithm history

    #     .. note::
    #         Because models are large, and multiple iterations of models need to
    #         be stored in memory, previous models are stored as `memmaps`,
    #         which allow for access of small segments of large files on disk,
    #         without reading the entire file. Memmaps are array like objects.

    #     .. note::
    #         Notation for s and y taken from Liu & Nocedal 1989
    #         iterate notation: sk = x_k+1 - x_k and yk = g_k+1 - gk

    #     :rtype s: np.memmap
    #     :return s: memory of the model differences `m_new - m_old`
    #     :rtype y: np.memmap
    #     :return y: memory of the gradient differences `g_new - g_old`
    #     """
    #     # Determine the iterates for model m and gradient g
    #     s_k = \
    #         self.load_vector("m_new").vector - self.load_vector("m_old").vector
    #     y_k = \
    #         self.load_vector("g_new").vector - self.load_vector("g_old").vector

    #     # Determine the shape of the memory map (length of model, length of mem)
    #     m = len(s_k)
    #     n = self.LBFGS_mem

    #     # Initial iteration, need to create the memory map
    #     if self._memory_used == 0:
    #         s = np.memmap(filename=self.path._s_file, mode="w+",
    #                       dtype="float32", shape=(m, n))
    #         y = np.memmap(filename=self.path._y_file, mode="w+",
    #                       dtype="float32", shape=(m, n))
    #         # Store the model and gradient differences in memmaps
    #         s[:, 0] = s_k
    #         y[:, 0] = y_k
    #         self._memory_used = 1
    #     # Subsequent iterations will append to memory maps
    #     else:
    #         s = np.memmap(filename=self.path._s_file, mode="r+",
    #                       dtype="float32", shape=(m, n))
    #         y = np.memmap(filename=self.path._y_file, mode="r+",
    #                       dtype="float32", shape=(m, n))
    #         # Shift all stored memory by one index to make room for latest mem
    #         s[:, 1:] = s[:, :-1]
    #         y[:, 1:] = y[:, :-1]
    #         # Store the latest model and gradient in first index
    #         s[:, 0] = s_k
    #         y[:, 0] = y_k

    #         # Keep track of the memory used
    #         if self._memory_used < self.LBFGS_mem:
    #             self._memory_used += 1

    #     return s, y

    def _apply_inverse_hessian(self, q, s=None, y=None):
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
        # If no memmaps are given as arguments, instantiate them
        if s is None or y is None:
            m = len(q)
            n = self.LBFGS_mem
            s = np.memmap(filename=self.path._s_file, mode="w+",
                          dtype="float32", shape=(m, n))
            y = np.memmap(filename=self.path._y_file, mode="w+",
                          dtype="float32", shape=(m, n))

        # First matrix product
        # Recursion step 2 from appendix A of Modrak & Tromp 2016
        kk = self._memory_used
        rh = np.zeros(kk)
        al = np.zeros(kk)
        for ii in range(kk):
            rh[ii] = 1 / np.dot(y[:, ii], s[:, ii])
            al[ii] = rh[ii] * np.dot(s[:, ii], q)
            q = q - al[ii] * y[:, ii]

        # Apply an optional preconditioner. Otherwise r==q
        # !!! MODEL - I need to rethink how this is going to work
        if self.preconditioner:
            r = self.precondition(q)
        else:
            # This is just a reference so no memory copy
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
