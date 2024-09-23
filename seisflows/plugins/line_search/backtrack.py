#!/usr/bin/env python3
"""
Backtracking line search class plugin to be used with an L-BFGS optimization

https://en.wikipedia.org/wiki/Backtracking_line_search
"""
from seisflows import logger
from seisflows.plugins.line_search.bracket import Bracket
from seisflows.tools.math import parabolic_backtrack


class Backtrack(Bracket):
    """
    [line_search.backtrack] Backtracking line search assumes the
    gradient is well scaled from the L-BFGS optimization algorithm, such
    that a unit step length (1) will provide a decrease in misfit. If
    misfit does not decrease, the backtracking step count follows a
    parabolic backtrack from 1 -> 0 in search of a decreased misfit. If the
    backtracked value becomes too small the backtracking line search defaults to
    a Bracketing line search.

    Variables Descriptions:
        x: list of step lenths from current line search
        f: correpsonding list of function values
        m: number of step lengths in current line search
        n: number of model updates in optimization problem
        gtg: dot product of gradient with itself
        gtp: dot product of gradient and search direction

    Status codes
        status == 1   : PASS, line search finished
        status == 0   : TRY/RETRY, attempt line search w/ new step length
        status == -1  : FAIL, line search exceeds internal criteria
    """
    def calculate_step_length(self):
        """
        Determines step length and search status. Defaults to 'Bracket'ing
        line search during the first evaluation (waiting for the L-BFGS to scale
        properly).

        .. note::
            Search history variable descriptions:
            x: list of step lenths from current line search
            f: correpsonding list of function values
            m: number of step lengths in current line search
            n: number of model updates in optimization problem
            gtg: dot product of gradient with itself
            gtp: dot product of gradient and search direction

        .. note:
            Available status returns are:
            'TRY': try/re-try the line search as conditions have not been met
            'PASS': line search was successful, you can terminate the search
            'FAIL': line search has failed for one or more reasons.

        :rtype: tuple (float, str)
        :return: (alpha==calculated step length,
            status==how to treat the next step count evaluation)
        """
        # Determine the current line search history
        x, f = self.get_search_history()
        first_iteration = bool(self.update_count == 1)
        first_step = bool(self.step_count == 1)
        
        # Check if this is the first iteration (by counting 0 step lens).
        # Quasi-Newton direction is not yet scaled properly if first iter. 
        # Instead of a bactracking line perform a bracketing line search
        if first_iteration:
            logger.info("first iteration, defaulting to bracketing line search")
            alpha, status = super().calculate_step_length()
        # Assumed well scaled search direction, attempt backtracking line search 
        # with unit step length
        else:
            # Initial unit step length
            if first_step:
                alpha = 1.
                logger.info(f"try: first step of iteration, attempt unit step")
                status = "TRY"
            # Pass if misfit is reduced
            elif f.min() < f[0]:
                alpha = x[f.argmin()]
                logger.info(f"pass: misfit decreased, line search successful")
                status = "PASS"
            # If misfit continually increases, decrease step length
            elif self.step_count < self.step_count_max:
                slope = self.gtp[-1] / self.gtg[-1]
                alpha = parabolic_backtrack(f0=f[0], g0=slope, x1=x[1],
                                            f1=f[1], b1=0.1, b2=0.5)
                logger.info(f"try: misfit increasing, parabloic backtrack "
                            f"to decrease step length")
                status = "TRY"
            # Failed because step_count_max exceeded
            else:
                logger.info(f"fail: backtracking line search has "
                            f"failed because the maximum allowable step counts "
                            f"({self.step_count_max}) has been exceeded")
                alpha = None
                status = "FAIL"
            # Decrement step count to set the final accepted value if pass 
            # because we will not be running another step
            if status == "PASS":
                self.step_count -= 1
                logger.info(f"final accepted step count == {self.step_count}")

        return alpha, status

