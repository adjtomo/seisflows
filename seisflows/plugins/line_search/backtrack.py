#!/usr/bin/env python3
"""
This is the subclass class for seisflows.plugins.line_search.backtrack
"""
from seisflows.plugins.line_search.bracket import Bracket
from seisflows.tools import msg
from seisflows.tools.math import parabolic_backtrack


class Backtrack(Bracket):
    """
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
    def __init__(self):
        """Inherits from seisflows.plugins.line_search.bracket.Bracket"""
        super().__init__()

    def initialize(self, step_len, func_val, gtg, gtp):
        """Inherits from seisflows.plugins.line_search.bracket.Bracket"""
        self.initialize(step_len, func_val, gtg, gtp)

    def update(self, step_len, func_val):
        """Inherits from seisflows.plugins.line_search.bracket.Bracket"""
        self.update(step_len, func_val)

    def clear_history(self):
        """Inherits from seisflows.plugins.line_search.bracket.Bracket"""
        self.clear_history()

    def reset(self):
        """Inherits from seisflows.plugins.line_search.bracket.Bracket"""
        self.reset()

    def search_history(self, sort=True):
        """Inherits from seisflows.plugins.line_search.bracket.Bracket"""
        return self.search_history(sort=sort)

    def calculate_step(self):
        """
        Determines step length and search status. Overwrites the Bracketing
        line search step calculation
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
            self.logger.debug(msg.sub(f"BACKTRACKING LINE SEARCH STEP"
                                      f"{self.step_count:0>2}"))
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
                import pdb;pdb.set_trace()
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

