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
        self.iteridx = []

    def __str__(self):
        """Simply print out all attributes, used for debugging mainly"""
        str_out = ""
        for key, val in vars(self).items():
            str_out += f"{key}: {val}\n"
        return str_out
    
    @property
    def update_count(self):
        """Independently keeps track of how many model updates we have made"""
        return len(self.iteridx)

    def check(self):
        """
        Make sure internal parameters are consistent with expectations. This is
        always expected to pass if the User/workflow uses the main functions 
        to manipulate the line search. But occasionally some manual intervention
        is required, so the check function makes sure things are okay.
        """
        try:
            assert len(self.iteridx) == len(self.gtg), (
                "number of iterations is inconsistent with intializations"
            )
            assert(len(self.gtg) == len(self.gtp)), (
                "line search GTP and GTG have become inconsistent"
            )
            assert(len(self.step_lens) == len(self.func_vals)), (
                "line search step lengths and func vals inconsistent"
            )
        except AssertionError as e: 
            logger.critical(e)

    def initialize_line_search(self, func_val, gtg, gtp):
        """
        Input the initial values of a line search, corresponding to the 
        initial model for this given iteration. Step length is assumed to be 0,
        because step length is calculated relative to the initial model.

        .. note::

            This is a one-time permanent change to the search history, but can 
            be rolled back manually by the User with `_restart_line_search` (to 
            start over from before `initialize_line_search`

        :type func_val: float
        :param func_val: value of the objective function evaluation
        :type step_len: float
        :param step_len: length of step for the trial model (alpha)
        :type gtg: float
        :param gtg: dot product of gradient with itself
        :type gtp: float
        :param gtp: dot product of gradient with search direction
        """
        # Append information regarding the current line search
        self.step_count = 1
        self.step_lens.append(0)
        self.func_vals.append(func_val)
        self.gtg.append(gtg)
        self.gtp.append(gtp)

        # `iteridx`` keeps track of the index of when each iteration started 
        # based on the step length. That way if we need to revert or restart a 
        # line search, we know exactly where to roll back to
        self.iteridx.append(len(self.step_lens) - 1)  # -1 for idx, not length

    def update_search_history(self, func_val, step_len):
        """
        After a misfit evaluation for a trial model, line search needs to know
        how large the step (alpha) was, and the corresponding objective function
        misfit evaluation value.

        .. note::

            This is a one-time permanent change to the search history, but can 
            be rolled back manually by the User with `_revert_line_search`, to 
            undo the update performed here

        :type func_val: float
        :param func_val: value of the objective function evaluation
        :type step_len: float
        :param step_len: length of step for the trial model (alpha)
        """
        # Quick check to make sure this is expected
        x, f = self.get_search_history()
        assert len(x) == self.step_count, \
            "update was not expected, step count does not match array length"
        
        self.func_vals.append(func_val)
        self.step_lens.append(step_len)

    def clear_search_history(self):
        """
        Clears internal line search history after an optimization restart. This
        is the 'nuclear' option, because it gets rid of ALL saved information
        about misfit evaluations from the previous iterations, i.e., this is 
        like starting from scratch.
        """
        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []
        self.step_count = 0
        self.iteridx = []

    def _restart_line_search(self):
        """
        Iff a line search has been initialized and something goes awry, it is
        often preferable to restart the line search. To do this, we roll back
        the number of steps we have taken during the line search, and undo
        variable changes that took place in 'initialize_search'. 
        
        .. note::
            In order to proceed with workflow after running this function you 
            will need to restart the workflow from `initialize_line_search`
        """
        logger.info("restarting line search for the current iteration")

        # Incase we somehow get doubles which should never be the case
        self.iteridx = np.unique(self.iteridx).tolist()
        self.step_count = 0

        # Clear out search history up to the last iteration
        idx = self.iteridx[-1] 
        self.step_lens = self.step_lens[:idx]
        self.func_vals = self.func_vals[:idx]

        # Roll back one-time initialized parameters to before initialization
        idx = self.update_count - 1  # roll back to the last time we updated
        self.iteridx = self.iteridx[:idx]
        self.gtg = self.gtg[:idx]
        self.gtp = self.gtp[:idx]

    def _revert_search_history(self):
        """
        Occasionally a line search will break mid-search but the User doesn't 
        want to `restart_line_search`, but rather just roll back to the previous
        step count. This function steps back the search history by 
        a number by one step. Call it multiple times to step back multiple steps
        Raises an error if there are no more steps to revert. 
        """
        if self.step_count == 0:
            logger.warning("cannot revert search history, already at step "
                           "count 0. Use `restart_search_history` to restart "
                           "to a pre line search state")
            return
        
        logger.info(f"reverting search history 1 step")

        # Step back the other quantities however many steps we have taken
        self.func_vals = self.func_vals[:-1]
        self.step_lens = self.step_lens[:-1]
        self.step_count -= 1

    def get_search_history(self):
        """
        Get the step lengths and misfit function evaluations for the current
        line search. Sort the list by the step lengths to get a coherent list,
        as e.g., in a bracketing line search the actual order of the list 
        does not match the magnitude of step lengths taken.

        .. note::

            In the original SeisFlows the values of the search history were
            non-intuitive and difficult to understand at a glance. The 
            original values are copied here as they relate to the mathematical
            formulations and so are still important:
        
            i = self.step_count  
            k = len(self.step_lens)
            x = np.array(self.step_lens[k - i - 1:k])
            f = np.array(self.func_vals[k - i - 1:k])
            j = count_zeros(self.step_lens) - 1  # update count

        :rtype: (np.array, np.array)
        :return: (step lengths of current line search, 
                  misfits of current line search)
        """
        if not self.iteridx:
            logger.warning("line search has not yet been initialized")
            return
        
        # Last iteration marks the index of the current line search
        idx = self.iteridx[-1]

        # Only return the last valid step lengths and misfit values that 
        # represent the current line search
        step_lens = np.array(self.step_lens[idx:])  # x
        func_vals = np.array(self.func_vals[idx:])  # f

        # Sort by the step lengths taken
        func_vals = func_vals[abs(step_lens).argsort()]
        step_lens = step_lens[abs(step_lens).argsort()]

        return step_lens, func_vals

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
        x, f = self.get_search_history()  
        # Some boolean checks to see where we're at in the inversion
        first_iteration = bool(self.update_count == 1)
        first_step = bool(self.step_count == 1)

        # For the first inversion and initial step, set alpha manually
        if first_iteration and first_step:  # == i00s00
            # Based on idea from Dennis and Schnabel
            alpha = self.gtg[-1] ** -1
            logger.info(f"try: first evaluation, guessing step length based on "
                        f"gradient value")
            status = "TRY"
        # For every iteration's initial step, set alpha manually
        elif first_step and (not first_iteration):
            # Based on the first equation in sec 3.5 of Nocedal and Wright 2ed
            idx = np.argmin(self.func_vals[:-1])
            alpha = self.step_lens[idx] * self.gtp[-2] / self.gtp[-1]
            logger.info(f"try: first step of iteration, "
                        f"setting scaled step length")
            status = "TRY"
        # If misfit is reduced and then increased, we've bracketed. Pass
        elif _check_bracket(x, f) and _good_enough(x, f):
            alpha = x[f.argmin()]
            logger.info(f"pass: bracket acceptable and step length "
                        f"reasonable. returning minimum line search misfit")
            status = "PASS"
        # If misfit is reduced but not close, set to quadratic fit
        elif _check_bracket(x, f):
            alpha = polynomial_fit(x, f)
            logger.info(f"try: bracket acceptable but step length unreasonable "
                        f"attempting to re-adjust step length")
            status = "TRY"
        # If misfit continues to step down, increase step length
        elif self.step_count < self.step_count_max and all(f <= f[0]):
            alpha = 1.618034 * x[-1]  # 1.618034 is the 'golden ratio'
            logger.info(f"try: misfit not bracketed, increasing step length "
                        f"using golden ratio")
            status = "TRY"
        # If misfit increases, reduce step length by backtracking
        elif self.step_count < self.step_count_max:
            slope = self.gtp[-1] / self.gtg[-1]
            alpha = parabolic_backtrack(f0=f[0], g0=slope, x1=x[1],
                                        f1=f[1], b1=0.1, b2=0.5)
            logger.info(f"try: misfit increasing, attempting "
                        f"to reduce step length using parabolic backtrack")
            status = "TRY"
        # step_count_max exceeded, fail
        else:
            logger.info(f"fail: bracketing line search has failed "
                        f"to reduce the misfit before exceeding "
                        f"`step_count_max`={self.step_count_max}")
            alpha = None
            status = "FAIL"
                
        # Decrement step count to set the final accepted value if pass because
        # we will not be running another step
        if status == "PASS":
            self.step_count -= 1
            logger.info(f"final accepted step count == {self.step_count}")

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




