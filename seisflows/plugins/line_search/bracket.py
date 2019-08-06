
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

        print '\tBracketing line search'
        print '\t\tStep lengths = {}'.format(x)
        print '\t\tMisfits = {}'.format(f)
        
        # For the first inversion and initial step, set alpha
        if step_count==0 and update_count==0:
            # based on idea from Dennis and Schnabel
            print "\t\tFirst inversion, initial trial step, continuing..."
            alpha = gtg[-1]**-1
            status = 0

        # For every i'th inversion's initial step, set alpha
        elif step_count==0:
            # based on the first equation in sec 3.5 of Nocedal and Wright 2ed
            print "\t\tInitial trial step, setting scaled step length"
            idx = np.argmin(self.func_vals[:-1])
            alpha = self.step_lens[idx] * gtp[-2]/gtp[-1]
            status = 0

        # If misfit is reduced and then increased, we've bracketed. Pass 
        elif _check_bracket(x,f) and _good_enough(x,f):
            print "\t\tBracket okay, step length reasonable, pass" 
            alpha = x[f.argmin()]
            status = 1

        # if misfit is reduced but not close, set to quadratic fit
        elif _check_bracket(x,f):
            print "\t\tBracket okay, step length unreasonable, manual step..."
            alpha = polyfit2(x,f)
            status = 0

        # if misfit continues to step down, increase step length
        elif step_count <= self.step_count_max and all(f <= f[0]):
            print "\t\tMisfit not bracketed, increasing step length..."
            alpha = 1.618034*x[-1]
            status = 0

        # if misfit increases, reduce step length
        elif step_count <= self.step_count_max:
            print "\t\tMisfit increasing, reducing step length..."
            slope = gtp[-1]/gtg[-1]
            alpha = backtrack2(f[0], slope, x[1], f[1], b1=0.1, b2=0.5)
            status = 0

        # step_count_max exceeded, fail
        else:
            print "\t\tBracketing failed, step_count_max exceeded"
            alpha = None
            status = -1

        # apply optional step length safeguard
        if alpha > self.step_len_max and step_count==0:
            print "\tInitial step length safegaurd, setting manual step length"
            alpha = 0.618034*self.step_len_max
            status = 0
            
        # stop because safeguard prevents us from going further
        elif alpha > self.step_len_max:
            print "step_len_max exceeded, manual set alpha"
            alpha = self.step_len_max
            status = 1

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
    :rtype: int
    :return: status of function as a bool
    """
    x, f = step_lens, func_vals
    imin, fmin = f.argmin(), f.min()
    if (fmin < f[0]) and any(f[imin:] > fmin):
        return 1
    else:
        return 0


def _good_enough(step_lens, func_vals, thresh=np.log10(1.2)):
    """
    Checks if step length is reasonably close to quadratic estimate
    
    polyfit2 requires 
    
    :type step_lens: np.array
    :param step_lens: an array of the step lengths taken during iteration
    :type func_vals: np.array
    :param func_vals: array of misfit values from eval func function
    :type thresh: numpy.float64
    :param thresh: threshold value for comparison against quadratic estimate
    :rtype: int
    :return: status of function as a bool
    """
    x, f = step_lens, func_vals
    if not _check_bracket(x,f):
        return 0
    x0 = polyfit2(x,f)
    if any(np.abs(np.log10(x[1:]/x0)) < thresh):
        return 1
    else:
        return 0




