#!/usr/bin/env python
"""
This is the subclass class for seisflows.plugins.line_search.backtrack
"""
from seisflows.plugins.line_search.bracket import Bracket
from seisflows.tools.math import backtrack2


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

    def calculate_step(self):
        """
        Determines step length and search status
        """
        # Determine the line search history
        x, f, gtg, gtp, step_count, update_count = self.search_history()
        
        # quasi-Newton direction is not yet scaled properly, so instead
        # of a bactracking line perform a bracketing line search
        if update_count == 0:
            alpha, status = super(Backtrack, self).calculate_step()
       
        # Assumed well scaled search direction, attempt backtracking line search 
        # with unit step length
        else:
            if self.verbose:
                print("\tBacktracking line search")
                print(f"\t\tStep lengths = {x}")
                print(f"\t\tMisfits = {f}")
            # Initial unit step length
            if step_count == 0:
                if self.verbose:
                    print("\t\tAttempting unit step length")
                alpha = min(1., self.step_len_max)
                status = 0
            # Pass if misfit is reduced
            elif _check_decrease(x, sf):
                if self.vervose:
                    print("\t\tMisfit decrease, pass")
                alpha = x[f.argmin()]
                status = 1
            # If misfit continually increases, decrease step length
            elif step_count <= self.step_count_max:
                if self.verbose:
                    print("\t\tMisfit increase, decreasing step length")
                slope = gtp[-1] / gtg[-1]
                alpha = backtrack2(f0=f[0], g0=slope, x1=x[1], f1=f[1], b1=0.1,
                                   b2=0.5)
                status = 0
            # Failed because step_count_max exceeded
            else:
                if self.verbose:
                    print("\t\tBacktracking failed, step_count_max exceeded")
                alpha = None
                status = -1

        return alpha, status


def _check_decrease(step_lens, func_vals, c=1.e-4):
    """
    Checks for sufficient decrease by comparing the current functional value
    with the smallest functional value in the list.

    !!! This function is defined outside the class because the Base class
    !!! doesn't require the function, but this subclass does
    """
    x, f = step_lens, func_vals
    if f.min() < f[0]:
        return 1
    else:
        return 0

