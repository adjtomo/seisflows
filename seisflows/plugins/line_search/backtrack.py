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
        # Determine the line search history
        x, f, gtg, gtp, step_count, update_count = self.get_search_history()

        # quasi-Newton direction is not yet scaled properly, so instead
        # of a bactracking line perform a bracketing line search
        if update_count == 0:
            alpha, status = super().calculate_step_length()
       
        # Assumed well scaled search direction, attempt backtracking line search 
        # with unit step length
        else:
            self._print_stats(x, f)

            # Initial unit step length
            if step_count == 0:
                alpha = min(1., self.step_len_max)
                logger.info(f"try: attempt unit step length w/ alpha={alpha}")
                status = "TRY"
            # Pass if misfit is reduced
            elif f.min() < f[0]:
                alpha = x[f.argmin()]
                logger.info(f"pass: misfit decreased, line search "
                            f"successful w/ alpha={alpha}")
                status = "PASS"
            # If misfit continually increases, decrease step length
            elif step_count < self.step_count_max:
                slope = gtp[-1] / gtg[-1]
                alpha = parabolic_backtrack(f0=f[0], g0=slope, x1=x[1],
                                            f1=f[1], b1=0.1, b2=0.5)
                logger.info(f"try: misfit increasing, attempting "
                            f"to decrease step length to alpha={alpha}")
                status = "TRY"
            # Failed because step_count_max exceeded
            else:
                logger.info(f"fail: backtracking line search has "
                            f"failed because the maximum allowable step counts "
                            f"({self.step_count_max}) has been exceeded")
                alpha = None
                status = "FAIL"

        return alpha, status

