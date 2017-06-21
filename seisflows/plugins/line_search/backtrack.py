
from seisflows.plugins.line_search import Bracket
from seisflows.tools.math import backtrack2

import numpy as np


class Backtrack(Bracket):
    """ Implements backtracking linesearch
    """

    def initial_step(self):
        alpha = 1.
        if alpha > self.step_len_max:
            alpha = self.step_len_max
        return alpha


    def update(self):
        """ Checks if termination conditions are satisfied, and if necessary
          determines next step length in line search

          status > 0  : finished
          status == 0 : not finished
          status < 0  : failed
        """
        x, f, m, n, gtg, gtp = self.current_vals()

        if n == 1:
            alpha, status = super(Backtrack, self).update()

        elif _check_decrease(x,f):
            alpha = x[f.argmin()]
            status = 1

        elif m <= self.step_count_max:
            slope = gtp[-1]/gtg[-1]
            alpha = backtrack2(f[0], slope, x[1], f[1], b1=0.1, b2=0.5)
            status = 0

        else:
            alpha = None
            status = -1

        return alpha, status



def _check_decrease(step_lens, func_vals, c=1.e-4):
    """ Checks for sufficient decrease
    """
    x,f = step_lens, func_vals
    if f.min() < f[0]:
        return 1
    else:
        return 0

