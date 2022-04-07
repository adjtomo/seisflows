#!/usr/bin/env python3
"""
This is the main class for seisflows.line_search.preconds.diagonal
This class provides the utilities for a diagonal preconditioner
"""
import os
import sys


class Diagonal(object):
    """
    User supplied diagonal preconditioner
    Rescales model parameters based on user supplied weights
    """
    def __init__(self):
        """
        Loads any required dependencies
        """
        PATH = sys.modules["seisflows_paths"]
        solver = sys.modules["seisflows_solver"]

        if "PRECOND" not in PATH:
            raise Exception

        if not os.path.exists(PATH.PRECOND):
            raise Exception

        self.path = PATH.PRECOND
        self.load = solver.load
        self.merge = solver.merge

    def __call__(self, q):
        """
        Applies preconditioner to given vector

        :type q: np.array
        :param q: search direction
        :rtype: np.array
        :return: preconditioned search direction
        """
        p = self.merge(self.load(self.path))
        return p * q


