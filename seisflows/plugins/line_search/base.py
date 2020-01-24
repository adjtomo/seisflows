#!/usr/bin/env python
"""
This is the Base class for seisflows.plugins.line_search

Line search is called on by the optimization procedure
"""
from os.path import abspath
from seisflows.tools.array import count_zeros

import numpy as np


class Base(object):
    """
    Abstract base class for line search

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
    def __init__(self, step_count_max=10, step_len_max=np.inf,
                 path=abspath("."), verbose=True):
        # Set maximum number of trial steps
        self.step_count_max = step_count_max

        # Optional maximum step length safeguard
        self.step_len_max = step_len_max

        # Prepare output log
        self.writer = Writer(path)

        # Prepare lists for line search history
        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []
        self.step_count = 0
        self.verbose = verbose

    def initialize(self, step_len, func_val, gtg, gtp):
        """
        Initialize a new search from step count 0

        :param step_len:
        :param func_val:
        :param gtg:
        :param gtp:
        :return:
        """
        self.step_count = 0
        self.step_lens += [step_len]
        self.func_vals += [func_val]
        self.gtg += [gtg]
        self.gtp += [gtp]

        self.writer(step_len, func_val)

        return self.calculate_step()

    def update(self, step_len, func_val):
        """
        Update search history
        :param step_len:
        :param func_val:
        :return:
        """
        self.step_count += 1
        self.step_lens += [step_len]
        self.func_vals += [func_val]

        self.writer(step_len, func_val)

        return self.calculate_step()

    def clear_history(self):
        """
        Clears line search history
        """
        self.func_vals = []
        self.step_lens = []
        self.gtg = []
        self.gtp = []

    def search_history(self, sort=True):
        """
        A convenience function, collects information needed to determine
        search status and calculate step length

        :type sort: bool
        :param sort: sort the search history by step length
        :rtype x: np.array
        :return x: list of step lenths from current line search
        :rtype f: np.array
        :return f: correpsonding list of function values
        :rtype gtg: list
        :return gtg: dot product dot product of gradient with itself
        :rtype gtp: list
        :return gtp: dot product of gradient and search direction
        :rtype i: int
        :return i: step_count
        :rtype j: int
        :return j: number of step lengths that are zero
        """
        i = self.step_count
        j = count_zeros(self.step_lens) - 1
        k = len(self.step_lens)
        x = np.array(self.step_lens[k - i - 1:k])
        f = np.array(self.func_vals[k - i - 1:k])

        # Sort by step length
        if sort:
            f = f[abs(x).argsort()]
            x = x[abs(x).argsort()]

        return x, f, self.gtg, self.gtp, i, j

    def calculate_step(self):
        """
        Determines step length and search status

        !!! Must be implemented by subclass !!!
        """
        raise NotImplementedError("Must be implemented by subclass")


class Writer(object):
    """
    Utility for writing one or more columns to text file.
    Used to write the line search history into a text file with a set format.
    """
    def __init__(self, path="./output.optim"):
        """
        Initiate the Writer class

        :type path: str
        :param path: path to the file that the writer will write to
        """
        self.iter = 0
        self.filename = abspath(path)
        self.write_header()

    def __call__(self, steplen=None, funcval=None):
        """
        When the function is called, it writes to self.filename

        :type steplen: float
        :param steplen: step length
        :type funcval: float
        :param funcval: misfit functional value for given step
        :return:
        """
        iter_ = "{iteration:10d}  {step_length:10.3e}  {function_value:10.3e}\n"
        step_ = "{space:s}  {step_length:10.3e}  {function_value:10.3e}\n"

        with open(self.filename, "a") as fileobj:
            # First iteration or step length of 0 means new iteration
            if self.iter == 0 or steplen == 0:
                self.iter += 1
                fileobj.write(iter_.format(iteration=self.iter,
                                           step_length=steplen,
                                           function_value=funcval)
                              )
            # Non-new iteration means trial step lengths, do not iterate
            else:
                fileobj.write(step_.format(space=" " * 10,
                                           step_length=steplen,
                                           function_value=funcval)
                              )

    def write_header(self):
        """
        Write the header of the text file
        Allows for extra headers to be written into file, call function will
        need to be changed
        """
        headers = ["ITER", "STEPLEN", "MISFIT"]

        with open(self.filename, "a") as fileobj:
            # Write the headers
            for header in headers:
                fileobj.write(f"{header:>10s}  ")
            fileobj.write('\n')
            # Write some separators
            for _ in range(len(headers)):
                separator = "=" * 10
                fileobj.write(f"{separator:>10s}  ")

            fileobj.write('\n')

    def newline(self):
        """
        Write a new line to the text file
        """
        with open(self.filename, "a") as fileobj:
                fileobj.write("\n")



