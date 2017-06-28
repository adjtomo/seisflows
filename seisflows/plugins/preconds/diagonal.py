
import sys
import numpy as np

from os.path import exists


class Diagonal(object):
    """ User supplied diagonal preconditioner

        Rescales model parameters based on user supplied weights
    """
    def __init__(self):
        """ Loads any required dependencies
        """
        PAR = sys.modules['seisflows_parameters']
        PATH = sys.modules['seisflows_paths']

        solver = sys.modules['seisflows_solver']

        if 'PRECOND' not in PATH:
            raise Exception

        if not exists(PATH.PRECOND):
            raise Exception

        self.path = PATH.PRECOND
        self.load = solver.load
        self.merge = solver.merge


    def __call__(self, q):
        """ Applies preconditioner to given vector
        """
        p = self.merge(self.load(self.path))
        return p*q


