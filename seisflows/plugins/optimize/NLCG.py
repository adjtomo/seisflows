#!/usr/bin/env python
"""
This is the main class for seisflows.line_search.optimize.NLCG
This class provides the core utilities for Non-linear conjugate gradient methodss

It is called in as a plugin for seisflows.optimize.NLCG
"""
import os
import numpy as np

from seisflows.tools import unix
from seisflows.tools.math import dot
from seisflows.tools.tools import loadtxt, savetxt, loadnpy, savenpy


class NLCG:
    """
    Nonlinear conjugate gradient method

    Optimization Variables:
        m: model
        f: objective function value
        g: gradient direction
        p: search direction

    Line Search Variables:
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
    def __init__(self, path='.', load=loadnpy, save=savenpy, thresh=1.,
                 maxiter=np.inf, precond=None, verbose=True):
        """
        Initialize the NLCG algorithm

        :type path: str
        :param path: path to the optization directory
        :type load: function
        :param load: function to use for loading optimization objects
        :type save: function
        :param save: function to use for saving optimization objects
        :type thresh: float
        :param thresh: a safeguard threshold based on the angle between the
            gradient direction and the new direction
        :type maxiter: int or np.inf
        :param maxiter: maximum number of iterations allowed before restarting
            the NLCG machinery
        :type precond: function
        :param precond: optional preconditioner function
        """
        self.path = path
        self.load = load
        self.save = save
        self.maxiter = maxiter
        self.thresh = thresh
        self.precond = precond
        self.verbose = verbose

        # Determine if NLCG has already been started
        try:
            self.iter = loadtxt(os.path.join(self.path, "NLCG", "iter"))
        except IOError:
            unix.mkdir(os.path.join(self.path, "NLCG"))
            self.iter = 0

    def __call__(self):
        """
        Returns NLCG search direction

        :rtype: tuple (np.array, int)
        :return: search direction, status of search
        """
        if self.verbose:
            print("\tComputing search direction w/ NLCG")

        self.iter += 1
        savetxt(os.path.join(self.path, "NLCG", "iter", self.iter))

        unix.cd(self.path)
        g_new = self.load("g_new")

        # If first iteration, search direction is the current gradient
        if self.iter == 1:
            return -g_new, 0
        # Force restart if the iterations have surpassed the maximum number iter
        elif self.iter > self.maxiter:
            if self.verbose:
                print("restarting NLCG... [periodic restart]")
            self.restart()
            return -g_new, 1

        # Compute search direction
        g_old = self.load("g_old")
        p_old = self.load("p_old")

        # Apply preconditioner and calculate beta
        if self.precond:
            beta = pollak_ribere(g_new, g_old, self.precond)
            p_new = -self.precond(g_new) + beta * p_old
        else:
            beta = pollak_ribere(g_new, g_old)
            p_new = -g_new + beta * p_old

        # Check restart conditions, return search direction and status
        if check_conjugacy(g_new, g_old) > self.thresh:
            if self.verbose:
                print("restarting NLCG... [loss of conjugacy]")
            self.restart()
            return -g_new, 1
        elif check_descent(p_new, g_new) > 0.:
            if self.verbose:
                print("restarting NLCG... [not a descent direction]")
            self.restart()
            return -g_new, 1
        else:
            return p_new, 0

    def restart(self):
        """
        Restarts NLCG algorithm
        """
        self.iter = 1
        savetxt(os.path.join(self.path, "NLCG", "iter", self.iter))


def fletcher_reeves(g_new, g_old, precond=lambda x: x):
    """
    One method for calculating beta in the NLCG Algorithm
    from Fletcher & Reeves, 1964

    :type g_new: np.array
    :param g_new: new search direction
    :type g_old: np.array
    :param g_old: old search direction
    :type precond: function
    :param precond: preconditioner, defaults to simple return
    :rtype: float
    :return: beta, the scale factor to apply to the old search direction to
        determine the new search direction
    """
    num = dot(precond(g_new), g_new)
    den = dot(g_old, g_old)
    beta = num / den

    return beta


def pollak_ribere(g_new, g_old, precond=lambda x: x):
    """
    One method for calculating beta in the NLCG Algorithm
    from Polak & Ribiere, 1969

    :type g_new: np.array
    :param g_new: new search direction
    :type g_old: np.array
    :param g_old: old search direction
    :type precond: function
    :param precond: preconditioner, defaults to simple return
    :rtype: float
    :return: beta, the scale factor to apply to the old search direction to
        determine the new search direction
    """
    num = dot(precond(g_new), g_new-g_old)
    den = dot(g_old, g_old)
    beta = num/den
    return beta


def check_conjugacy(g_new, g_old):
    """
    Check for conjugacy between two vectors

    :type g_new: np.array
    :param g_new: new search direction
    :type g_old: np.array
    :param g_old: old search direction
    :rtype: float
    :return: an element that proves conjugacy
    """
    return abs(dot(g_new, g_old) / dot(g_new, g_new))


def check_descent(p_new, g_new):
    """
    Ensure that the search direction is descending

    :type p_new: np.array
    :param p_new: search direction
    :type g_new: np.array
    :param g_new: gradient direction
    :rtype: float
    :return: the angle between search direction and gradient direction, should
        be negative to ensure descent
    """
    return dot(p_new, g_new) / dot(g_new, g_new)




