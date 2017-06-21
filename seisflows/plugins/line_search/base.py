

from seisflows.tools.array import count_zeros

import numpy as np


class Base(object):
    """ Abstract base class for line search

      Variables
          x - list of step lenths from current line search
          f - correpsonding list of function values
          m - how many step lengths in current line search?
          n - how many model updates in optimization problem?
          gtg - dot product of gradient with itself                    
          gtp - dot product of gradient and search direction

      Status codes
          status > 0  : finished
          status == 0 : not finished
          status < 0  : failed
    """
    def __init__(self,
                step_count_max=10,
                step_len_max=np.inf):

        # maximum number of trial steps
        self.step_count_max = step_count_max

        # optional maximum step length safeguard
        self.step_len_max = step_len_max

        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []


    def clear_history(self):
        """ Clears line search history
        """
        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []


    def current_vals(self, sort=True):
        """ Collects information about current line search
        """
        i = self.step_count
        j = count_zeros(self.step_lens)
        k = len(self.step_lens)
        x = np.array(self.step_lens[k-i-1:k])
        f = np.array(self.func_vals[k-i-1:k])
        if sort:
            f = f[abs(x).argsort()]
            x = x[abs(x).argsort()]
        return x, f, i, j, self.gtg, self.gtp


    def initial_step(self):
        # must be implemented by subclass
        raise NotImplementedError


    def update(self):
        # must be implemented by subclass
        raise NotImplementedError




