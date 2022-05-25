#!/usr/bin/env python3
"""
This is the subclass class for seisflows.plugins.line_search.bracket
"""
import logging
import numpy as np

from seisflows.tools import msg
from seisflows.plugins.line_search.base import Base
from seisflows.tools.math import parabolic_backtrack, polynomial_fit


class Bracket(Base):
    """
    Implements bracketing line search, which attempts to find a step length
    corresponding to misfit reduction, and a misfit corresponding to misfit
    increase, so that the optimization procedure can scale the step length
    in future iterations.

    Variables Descriptions:
        x: list of step lenths from current line search
        f: correpsonding list of function values
        gtg: dot product of gradient with itself
        gtp: dot product of gradient and search direction

    Status codes
        status > 0  : finished
        status == 0 : not finished
        status < 0  : failed
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self, **kwargs):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        super().__init__(**kwargs)

    def calculate_step(self):
        """
        Determines step length (alpha) and search status (status)
        """
        # Determine the line search history
        x, f, gtg, gtp, step_count, update_count = self.search_history()

        # Print out the current line search parameters for convenience
        self.logger.debug(msg.sub("EVALUATE BRACKETING LINE SEARCH"))
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
        elif self._check_bracket(x, f) and self._good_enough(x,f):
            alpha = x[f.argmin()]
            self.logger.info(f"bracket okay, step length reasonable, pass")
            status = 1
        # If misfit is reduced but not close, set to quadratic fit
        elif self._check_bracket(x, f):
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
                self.logger.info(f"step_len_max={self.step_len_max} exceeded, "
                                 f"manual set alpha")
                alpha = self.step_len_max
                status = 1

        return alpha, status

    @staticmethod
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

    def _good_enough(self, step_lens, func_vals, thresh=np.log10(1.2)):
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
        if not self._check_bracket(x, f):
            okay = False
        else:
            x0 = polynomial_fit(x, f)
            if any(np.abs(np.log10(x[1:] / x0)) < thresh):
                okay = True
            else:
                okay = False
        return okay




