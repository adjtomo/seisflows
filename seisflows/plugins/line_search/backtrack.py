#!/usr/bin/env python3
"""
This is the subclass class for seisflows.plugins.line_search.backtrack
"""
import logging

from seisflows.tools import msg
from seisflows.plugins.line_search.bracket import Bracket
from seisflows.tools.math import parabolic_backtrack


class Backtrack(Bracket):
    """
    Overwrites seisflows.plugins.line_search.Bracket

    Implements backtracking linesearch. A backtracking line search is used
    for L-BFGS optimization, where a unit step length is attempted, if this
    does not satisfy the misfit reduction criteria, the step length is
    `backtracked` to a smaller value. If the backtracked value becomes too small
    the backtracking line search defaults to a Bracketing line search.

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
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)


    def calculate_step(self):
        """
        Determines step length and search status. Overloads the bracketing
        line search
        """
        # Determine the line search history
        x, f, gtg, gtp, step_count, update_count = self.search_history()
        
        # quasi-Newton direction is not yet scaled properly, so instead
        # of a bactracking line perform a bracketing line search
        if update_count == 0:
            alpha, status = super().calculate_step()
       
        # Assumed well scaled search direction, attempt backtracking line search 
        # with unit step length
        else:
            self.logger.info(msg.sub("EVALUATE BACKTRACKING LINE SEARCH"))
            x_str = ", ".join([f"{_:.2E}" for _ in x])
            f_str = ", ".join([f"{_:.2E}" for _ in f])
            self.logger.debug(f"step length(s) = {x_str}")
            self.logger.debug(f"misfit val(s)  = {f_str}")

            # Initial unit step length
            if step_count == 0:
                self.logger.info("attempting unit step length")
                alpha = min(1., self.step_len_max)
                status = 0
            # Pass if misfit is reduced
            elif f.min() < f[0]:
                self.logger.info("misfit decrease, pass")
                alpha = x[f.argmin()]
                status = 1
            # If misfit continually increases, decrease step length
            elif step_count <= self.step_count_max:
                self.logger.info("misfit increase, decreasing step length")
                slope = gtp[-1] / gtg[-1]
                alpha = parabolic_backtrack(f0=f[0], g0=slope, x1=x[1],
                                            f1=f[1], b1=0.1, b2=0.5)
                status = 0
            # Failed because step_count_max exceeded
            else:
                self.logger.info("backtracking failed, step_count_max exceeded")
                alpha = None
                status = -1

        return alpha, status

