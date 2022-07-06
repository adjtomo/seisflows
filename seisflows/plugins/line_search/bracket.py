#!/usr/bin/env python3
"""
This is the Bracketing line search class for seisflows

Line search is called on by the optimization procedure and should not really
have any agency (i.e. it should not be able to iterate its own step count etc.,
this should be completely left to the optimization algorithm to keep everything
in one place)
"""
import logging
import numpy as np

from seisflows.core import Base
from seisflows.tools import msg
from seisflows.tools.array import count_zeros
from seisflows.tools.math import parabolic_backtrack, polynomial_fit


class Bracket(Base):
    """
    Abstract base class for line search

    Variables Descriptions:
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
    def __init__(self, step_count_max, step_len_max):
        """
        Initiate the line search machinery

        :type step_count_max: int
        :param step_count_max: maximum number of step counts before changing
            line search behavior. set by PAR.STEP_COUNT_MAX
        :type step_len_max: int
        :param step_len_max: maximum length of the step, defaults to infinity,
            that is unbounded step length. set by PAR.STEP_LEN_MAX
        """
        super().__init__()

        self.step_count_max = step_count_max
        self.step_len_max = step_len_max
        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []
        self.step_count = 0

    def initialize(self, step_len, func_val, gtg, gtp):
        """
        Initialize a new search from step count 0 and calculate the step
        direction and length

        :type iter: int
        :param iter: current iteration defined by OPTIMIZE.iter
        :type step_len: float
        :param step_len: initial step length determined by optimization
        :type func_val: float
        :param func_val: current evaluation of the objective function
        :type gtg: float
        :param gtg: dot product of the gradient with itself
        :type gtp: float
        :param gtp: dot product of gradient `g` with search direction `p`
        :rtype alpha: float
        :return alpha: the calculated trial step length
        :rtype status: int
        :return status: current status of the line search
        """
        self.step_count = 0
        self.step_lens.append(step_len)
        self.func_vals.append(func_val)
        self.gtg.append(gtg)
        self.gtp.append(gtp)

    def update(self, step_len, func_val):
        """
        Update search history by appending internal attributes, writing the
        current list of step lengths and function evaluations, and calculating a
        new step length

        :type iter: int
        :param iter: current iteration defined by OPTIMIZE.iter
        :type step_len: float
        :param step_len: step length determined by optimization
        :type func_val: float
        :param func_val: current evaluation of the objective function
        :rtype alpha: float
        :return alpha: the calculated rial step length
        :rtype status: int
        :return status: current status of the line search
        """
        self.step_lens.append(step_len)
        self.func_vals.append(func_val)

    def clear_history(self):
        """
        Clears internal line search history for a new line search attempt
        """
        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []
        self.step_count = 0

    def reset(self):
        """
        If a line search fails mid-search, and the User wants to resume from 
        the line search function. Initialize will be called again. This function
        undos the progress made by the previous line search so that a new line
        search can be called without problem.

        output.optim needs to have its lines cleared manually
        """
        # First step treated differently
        if len(self.step_lens) <= 1:
            self.clear_history()
        else:
            # Wind back dot products by one
            self.gtg = self.gtg[:-1]
            self.gtp = self.gtp[:-1]
            
            # Move step lens and function evaluations by number of step count
            original_idx = -1 * self.step_count - 1
            self.step_lens = self.step_lens[:original_idx]
            self.func_vals = self.func_vals[:original_idx]

    def search_history(self, sort=True):
        """
        A convenience function, collects information based on the current
        evaluation of the line search, needed to determine search status and 
        calculate step length. From the full collection of the search history,
        only returns values relevant to the current line search.

        :type sort: bool
        :param sort: sort the search history by step length
        :rtype x: np.array
        :return x: list of step lenths from current line search
        :rtype f: np.array
        :return f: correpsonding list of function values
        :rtype gtg: list
        :return gtg: dot product dot product of gradient with itself
        :rtype gtp: list
        :return gtp: dot product of gradient and search direction
        :rtype i: int
        :return i: step_count
        :rtype j: int
        :return j: number of iterations corresponding to 0 step length
        """
        i = self.step_count
        j = count_zeros(self.step_lens) - 1
        k = len(self.step_lens)
        x = np.array(self.step_lens[k - i - 1:k])
        f = np.array(self.func_vals[k - i - 1:k])

        # Sort by step length
        if sort:
            f = f[abs(x).argsort()]
            x = x[abs(x).argsort()]

        return x, f, self.gtg, self.gtp, i, j

    def calculate_step(self):
        """
        Determines step length (alpha) and search status (status)
        """
        # Determine the line search history
        x, f, gtg, gtp, step_count, update_count = self.search_history()

        # Print out the current line search parameters for convenience
        self.logger.debug(msg.sub(f"BRACKETING LINE SEARCH STEP "
                                  f"{self.step_count:0>2}"))
        x_str = ", ".join([f"{_:.2E}" for _ in x])
        f_str = ", ".join([f"{_:.2E}" for _ in f])
        self.logger.debug(f"step length(s) = {x_str}")
        self.logger.debug(f"misfit val(s)  = {f_str}")

        # For the first inversion and initial step, set alpha manually
        if step_count == 0 and update_count == 0:
            # Based on idea from Dennis and Schnabel
            alpha = gtg[-1] ** -1
            self.logger.info(f"first iteration, guessing trial step")
            status = 0
        # For every i'th inversions initial step, set alpha manually
        elif step_count == 0:
            # Based on the first equation in sec 3.5 of Nocedal and Wright 2ed
            idx = np.argmin(self.func_vals[:-1])
            alpha = self.step_lens[idx] * gtp[-2] / gtp[-1]
            self.logger.info(f"first step, setting scaled step length")
            status = 0
        # If misfit is reduced and then increased, we've bracketed. Pass
        elif _check_bracket(x, f) and _good_enough(x, f):
            alpha = x[f.argmin()]
            self.logger.info(f"bracket okay, step length reasonable, pass")
            status = 1
        # If misfit is reduced but not close, set to quadratic fit
        elif _check_bracket(x, f):
            alpha = polynomial_fit(x, f)
            self.logger.info(f"bracket okay, step length unreasonable, "
                             f"manual step")
            status = 0
        # If misfit continues to step down, increase step length
        elif step_count <= self.step_count_max and all(f <= f[0]):
            alpha = 1.618034 * x[-1]  # 1.618034 is the 'golden ratio'
            self.logger.info(f"misfit not bracketed, increasing step length")
            status = 0
        # If misfit increases, reduce step length by backtracking
        elif step_count <= self.step_count_max:
            slope = gtp[-1] / gtg[-1]
            alpha = parabolic_backtrack(f0=f[0], g0=slope, x1=x[1],
                                        f1=f[1], b1=0.1, b2=0.5)
            self.logger.info(f"misfit increasing, reducing step length to")
            status = 0
        # step_count_max exceeded, fail
        else:
            self.logger.info(f"bracketing failed, "
                             f"step_count_max={self.step_count_max} exceeded")
            alpha = None
            status = -1

        # Apply optional step length safeguard
        if alpha is not None:
            if alpha > self.step_len_max and step_count == 0:
                alpha = 0.618034 * self.step_len_max
                self.logger.info(f"initial step length safegaurd, setting "
                                 f"manual step length")
                status = 0
            # Stop because safeguard prevents us from going further
            elif alpha > self.step_len_max:
                self.logger.info(f"step_len_max={self.step_len_max:.2f} "
                                 f"exceeded, manual set alpha")
                alpha = self.step_len_max
                status = 1

        return alpha, status


def _check_bracket(step_lens, func_vals):
    """
    Checks if minimum has been bracketed

    Looks at the minimum of the misfit values calculated through eval func
    to see if the misfit has been reduced w.r.t the initial misfit

    :type step_lens: numpy.array
    :param step_lens: an array of the step lengths taken during iteration
    :type func_vals: numpy.array
    :param func_vals: array of misfit values from eval func function
    :rtype: bool
    :return: status of function as a bool
    """
    x, f = step_lens, func_vals
    imin, fmin = f.argmin(), f.min()
    if (fmin < f[0]) and any(f[imin:] > fmin):
        okay = True
    else:
        okay = False
    return okay


def _good_enough(step_lens, func_vals, thresh=np.log10(1.2)):
    """
    Checks if step length is reasonably close to quadratic estimate

    :type step_lens: np.array
    :param step_lens: an array of the step lengths taken during iteration
    :type func_vals: np.array
    :param func_vals: array of misfit values from eval func function
    :type thresh: numpy.float64
    :param thresh: threshold value for comparison against quadratic estimate
    :rtype: bool
    :return: status of function as a bool
    """
    x, f = step_lens, func_vals
    if not _check_bracket(x, f):
        okay = False
    else:
        x0 = polynomial_fit(x, f)
        if any(np.abs(np.log10(x[1:] / x0)) < thresh):
            okay = True
        else:
            okay = False
    return okay




