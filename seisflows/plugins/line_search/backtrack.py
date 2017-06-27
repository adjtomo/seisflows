
from seisflows.plugins.line_search import Bracket
from seisflows.tools.math import backtrack2

import numpy as np


class Backtrack(Bracket):
    """ Implements backtracking linesearch

      Variables
          x - list of step lenths from current line search
          f - correpsonding list of function values
          gtg - dot product of gradient with itself                    
          gtp - dot product of gradient and search direction

      Status codes
          status > 0  : finished
          status == 0 : not finished
          status < 0  : failed
    """

    def calculate_step(self):
        """ Determines step length and search status
        """
        x, f, gtg, gtp, step_count, update_count = self.search_history()


        if update_count==0:
            # quasi-Newton direction is not yet scaled properly, so instead
            # of a bactracking line perform a bracketing line search
            alpha, status = super(Backtrack, self).calculate_step()

        elif step_count==0:
            # our choice of a unit step length here assumes a well-scaled
            # search direction
            alpha = min(1., self.step_len_max)
            status = 0

        elif _check_decrease(x,f):
            alpha = x[f.argmin()]
            status = 1

        elif step_count <= self.step_count_max:
            # we need a smaller step length
            slope = gtp[-1]/gtg[-1]
            alpha = backtrack2(f[0], slope, x[1], f[1], b1=0.1, b2=0.5)
            status = 0

        else:
            # failed because step_count_max exceeded
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

