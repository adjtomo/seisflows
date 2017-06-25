
from seisflows.plugins.line_search import Base
from seisflows.tools.math import backtrack2, polyfit2

import numpy as np



class Bracket(Base):
    """ Implements bracketing line search

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

        if step_count==0 and update_count==0:
            # based on idea from Dennis and Schnabel
            alpha = gtg[-1]**-1
            status = 0

        elif step_count==0:
            # based on the first equation in sec 3.5 of Nocedal and Wright 2ed
            idx = np.argmin(self.func_vals[:-1])
            alpha = self.step_lens[idx] * gtp[-2]/gtp[-1]
            status = 0

        elif _check_bracket(x,f) and _good_enough(x,f):
            alpha = x[f.argmin()]
            status = 1

        elif _check_bracket(x,f):
            alpha = polyfit2(x,f)
            status = 0

        elif step_count <= self.step_count_max and all(f <= f[0]):
            # we need a larger step length
            alpha = 1.618034*x[-1]
            status = 0

        elif step_count <= self.step_count_max:
            # we need a smaller step length
            slope = gtp[-1]/gtg[-1]
            alpha = backtrack2(f[0], slope, x[1], f[1], b1=0.1, b2=0.5)
            status = 0

        else:
            # failed because step_count_max exceeded
            alpha = None
            status = -1

        # apply optional step length safeguard
        if alpha > self.step_len_max and \
           step_count==0:
            alpha = 0.618034*self.step_len_max
            status = 0

        elif alpha > self.step_len_max:
            # stop because safeguard prevents us from going further
            alpha = self.step_len_max
            status = 1

        return alpha, status


def _check_bracket(step_lens, func_vals):
    """ Checks if minimum has been bracketed
    """
    x, f = step_lens, func_vals
    imin, fmin = f.argmin(), f.min()
    if (fmin < f[0]) and any(f[imin:] > fmin):
        return 1
    else:
        return 0


def _good_enough(step_lens, func_vals, thresh=np.log10(1.2)):
    """ Checks if step length is reasonably close to quadratic estimate
    """
    x, f = step_lens, func_vals
    if not _check_bracket(x,f):
        return 0
    x0 = polyfit2(x,f)
    if any(np.abs(np.log10(x[1:]/x0)) < thresh):
        return 1
    else:
        return 0




