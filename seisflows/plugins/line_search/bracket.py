#!/usr/bin/env python3
"""
A bracketing line search (a.k.a direct line search) attempts to find an
appropriate step length by identifying two points between which the minimum
misfit lies. Contains some functionality for saving line search history to disk
so that a line search may be resumed in the case of a failure/reset.

https://en.wikipedia.org/wiki/Line_search

.. note::
    Line search is called on by the optimization procedure and should not
    really have any agency (i.e. it should not be able to iterate its own step
    count etc., this should be completely left to the optimization algorithm to
    keep everything in one place)
"""
import os
import numpy as np

from seisflows import logger
from seisflows.tools.array import count_zeros
from seisflows.tools.math import parabolic_backtrack, polynomial_fit


class Bracket:
    """
    [line_search.bracket] The bracketing line search identifies two points
    between which the minimum misfit lies between.

    :type step_count_max: int
    :param step_count_max: maximum number of step counts before changing
        line search behavior. set by PAR.STEP_COUNT_MAX
    :type step_len_max: int
    :param step_len_max: maximum length of the step, defaults to infinity,
        that is unbounded step length. set by PAR.STEP_LEN_MAX
    """
    def __init__(self, step_count_max, step_len_max, path=None):
        """
        Instantiate max criteria for line search
        """
        self.step_count_max = step_count_max
        self.step_len_max = step_len_max

        if path is None:
            self.path = os.path.join(os.getcwd(), "line_search")
        else:
            self.path = path

        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []
        self.step_count = 0

    def update_search_history(self, func_val, step_len, gtg=None, gtp=None):
        """
        Update the internal list of search history attributes. Lists like
        `func_vals` get appended to, while values like step_count are
        overwritten. Allowed to increment func_val and step_len by themselves
        """
        self.func_vals.append(func_val)
        self.step_lens.append(step_len)
        if gtg:
            self.gtg.append(gtg)
        if gtp:
            self.gtp.append(gtp)

    def clear_search_history(self):
        """
        Clears internal line search history for a new line search attempt
        """
        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []
        self.step_count = 0

    def check_search_history(self):
        """
        Since the line search is just a wrapper for list of numbers, check that
        search history hasn't been muddled up by ensuring that internal lists
        are the correct length for the given evaluation
        """
        assert(len(self.gtg) == len(self.gtp)), f"too many entries for 'gtg'"
        assert(len(self.func_vals) == len(self.step_lens)), \
            f"number of function evaluations does not match step lengths"
        if self.func_vals:
            assert(self.step_count + 1 == len(self.func_vals)), \
                f"current step count doesn't match the number of function evals"

    def get_search_history(self, sort=True):
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
        :return j: number of iterations corresponding to 0 step length,
            i.e., the update count
        """
        i = self.step_count
        k = len(self.step_lens)
        x = np.array(self.step_lens[k - i - 1:k])
        f = np.array(self.func_vals[k - i - 1:k])
        j = count_zeros(self.step_lens) - 1  # update count

        # Sort by step length
        if sort:
            f = f[abs(x).argsort()]
            x = x[abs(x).argsort()]

        return x, f, self.gtg, self.gtp, i, j

    def _print_stats(self, x, f):
        """Print out misfit values and step lengths to the logger"""
        x_str = ", ".join([f"{_:.2E}" for _ in x])
        f_str = ", ".join([f"{_:.2E}" for _ in f])
        logger.debug(f"step length(s) = {x_str}")
        logger.debug(f"misfit val(s)  = {f_str}")

    def calculate_step_length(self):
        """
        Determines step length (alpha) and search status (status) using a
        bracketing line search. Evaluates Wolfe conditions to determine if
        a step length is acceptable.

        .. note:
            Available status returns are:
            'TRY': try/re-try the line search as conditions have not been met
            'PASS': line search was successful, you can terminate the search
            'FAIL': line search has failed for one or more reasons.

        :rtype: tuple (float, str)
        :return: (alpha==calculated step length,
            status==how to treat the next step count evaluation)
        """
        # Determine the line search history
        x, f, gtg, gtp, step_count, update_count = self.get_search_history()
        self._print_stats(x, f)

        # For the first inversion and initial step, set alpha manually
        if step_count == 0 and update_count == 0:
            # Based on idea from Dennis and Schnabel
            alpha = gtg[-1] ** -1
            logger.info(f"try: first evaluation, attempt guess step length, "
                        f"alpha={alpha:.2E}")
            status = "TRY"
        # For every iteration's initial step, set alpha manually
        elif step_count == 0:
            # Based on the first equation in sec 3.5 of Nocedal and Wright 2ed
            idx = np.argmin(self.func_vals[:-1])
            alpha = self.step_lens[idx] * gtp[-2] / gtp[-1]
            logger.info(f"try: first step count of iteration, "
                        f"setting scaled step length, alpha={alpha:.2E}")
            status = "TRY"
        # If misfit is reduced and then increased, we've bracketed. Pass
        elif _check_bracket(x, f) and _good_enough(x, f):
            alpha = x[f.argmin()]
            logger.info(f"pass: bracket acceptable and step length "
                        f"reasonable. returning minimum line search misfit.")
            status = "PASS"
        # If misfit is reduced but not close, set to quadratic fit
        elif _check_bracket(x, f):
            alpha = polynomial_fit(x, f)
            logger.info(f"try: bracket acceptable but step length unreasonable "
                        f"attempting to re-adjust step length "
                        f"alpha={alpha:.2E}")
            status = "TRY"
        # If misfit continues to step down, increase step length
        elif step_count < self.step_count_max and all(f <= f[0]):
            alpha = 1.618034 * x[-1]  # 1.618034 is the 'golden ratio'
            logger.info(f"try: misfit not bracketed, increasing step length "
                        f"using golden ratio, alpha={alpha:.2E}")
            status = "TRY"
        # If misfit increases, reduce step length by backtracking
        elif step_count < self.step_count_max:
            slope = gtp[-1] / gtg[-1]
            alpha = parabolic_backtrack(f0=f[0], g0=slope, x1=x[1],
                                        f1=f[1], b1=0.1, b2=0.5)
            logger.info(f"try: misfit increasing, attempting "
                        f"to reduce step length using parabloic backtrack, "
                        f"alpha={alpha:.2E}")
            status = "TRY"
        # step_count_max exceeded, fail
        else:
            logger.info(f"fail: bracketing line search has failed "
                        f"to reduce the misfit before exceeding "
                        f"`step_count_max`={self.step_count_max}")
            alpha = None
            status = "FAIL"

        # Apply optional step length safeguard
        if alpha is not None:
            if alpha > self.step_len_max and step_count == 0:
                alpha = 0.618034 * self.step_len_max
                logger.info(f"try: applying initial step length "
                            f"safegaurd as alpha has exceeded maximum step "
                            f"length, alpha_new={alpha:.2E}")
                status = "TRY"
            # Stop because safeguard prevents us from going further
            # TODO Why is this passing? should status not be carried over?
            elif alpha > self.step_len_max:
                alpha = self.step_len_max
                logger.info(f"try: applying initial step length "
                            f"safegaurd as alpha has exceeded maximum step "
                            f"length, alpha_new={alpha:.2E}")
                status = "PASS"

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




