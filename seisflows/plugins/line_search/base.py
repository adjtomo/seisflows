

from os.path import abspath
from seisflows.tools.array import count_zeros
from seisflows.tools.seismic import StepWriter

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
                step_len_max=np.inf,
                path=abspath('.')):

        # maximum number of trial steps
        self.step_count_max = step_count_max

        # optional maximum step length safeguard
        self.step_len_max = step_len_max

        # prepare output log
        self.writer = StepWriter(path)

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


    def search_history(self, sort=True):
        """ A convenience function, collects information needed to determine
          search status and calculate step length
        """
        i = self.step_count
        j = count_zeros(self.step_lens)-1
        k = len(self.step_lens)
        x = np.array(self.step_lens[k-i-1:k])
        f = np.array(self.func_vals[k-i-1:k])
        if sort:
            f = f[abs(x).argsort()]
            x = x[abs(x).argsort()]
        return x, f, self.gtg, self.gtp, i, j


    def initialize(self, step_len, func_val, gtg, gtp):

        # update search history
        self.step_count = 0
        self.step_lens += [step_len]
        self.func_vals += [func_val]
        self.gtg += [gtg]
        self.gtp += [gtp]

        self.writer(step_len, func_val)

        return self.calculate_step()


    def update(self, step_len, func_val):

        # update search history
        self.step_count += 1
        self.step_lens += [step_len]
        self.func_vals += [func_val]

        self.writer(step_len, func_val)

        return self.calculate_step()


    def calculate_step(self):
        raise NotImplementedError('Must be implemented by subclass')




