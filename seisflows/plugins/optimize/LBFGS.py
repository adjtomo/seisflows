#!/usr/bin/env python
"""
This is the main class for seisflows.line_search.optimize.LBFGS
This class provides the core utilities for the Seisflows optimization schema.

It is called in as a plugin for seisflows.optimize.LBFGS
"""
import numpy as np

from seisflows.tools import unix
from seisflows.tools.tools import exists, loadnpy, savenpy
from seisflows.tools.math import angle


class LBFGS(object):
    """
    Limited-memory BFGS algorithm

    Includes optional safeguards: periodic restarting and descent conditions.

    To conserve memory, most vectors are read from disk rather than passed
    from a calling routine.

    L-BFGS Variables:
        s: memory of model differences
        y: memory of gradient differences

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

    def __init__(self, path=".", load=loadnpy, save=savenpy, memory=5,
                 thresh=0., maxiter=np.inf, precond=None, verbose=True):
        """
        Initialize the LBFGS algorithm

        :type path: str
        :param path: path to the optization directory
        :type load: function
        :param load: function to use for loading optimization objects
        :type save: function
        :param save: function to use for saving optimization objects
        :type memory: int
        :param memory: the number of previous iterations saved locally for
            generating a step direction and step length
        :type thresh: float
        :param thresh: a safeguard threshold based on the angle between the
            gradient direction and the new direction
        :type maxiter: int or np.inf
        :param maxiter: maximum number of iterations allowed before restarting
            the L-BFGS machinery
        :type precond: function
        :param precond: optional preconditioner function
        """
        # Create the LBFGS directory in the path
        assert exists(path)
        unix.cd(path)
        unix.mkdir("LBFGS")

        # Assign values
        self.path = path
        self.load = load
        self.save = save
        self.thresh = thresh
        self.maxiter = maxiter
        self.precond = precond
        self.memory = memory
        self.verbose = verbose

        self.iter = 0
        self.memory_used = 0

    def __call__(self):
        """
        Returns L-BFGS search direction

        :rtype: tuple (np.array, int)
        :return: search direction, status of search
        """
        if self.verbose:
            print("\tComputing search direction w/ L-BFGS")

        self.iter += 1
        unix.cd(self.path)

        # Load the current gradient direction, which is the L-BFGS search
        # direction if this is the first iteration
        g = self.load("g_new")
        if self.iter == 1:
            return -g, 0
        elif self.iter > self.maxiter:
            if self.verbose:
                print("\trestarting LBFGS... [periodic restart]")
            self.restart()
            return -g, 1

        # Update the search direction, apply the inverse Hessian
        # q becomes the new search direction g
        s, y = self.update()
        q = self.apply(g, s, y)

        # Determine if the new search direction is appropriate
        status = self.check_status(g, q)
        if status != 0:
            self.restart()
            return -g, status
        else:
            return -q, status

    def update(self):
        """
        Updates L-BFGS algorithm history

        Note:
            Because models are large, and multiple iterations of models need to
            be stored in memory, previous models are stored as `memmaps`,
            which allow for access of small segments of large files on disk,
            without reading the entire file. Memmaps are array like objects.

        Notation for s and y taken from Liu & Nocedal 1989
            iterate notation: sk = x_k+1 - x_k and yk = g_k+1 - gk

        :rtype s: np.memmap
        :return s: memory of the model differences `m_new - m_old`
        :rtype y: np.memmap
        :return y: memory of the gradient differences `g_new - g_old`
        """
        unix.cd(self.path)

        # Determine the iterates for model m and gradient g
        s_k = self.load("m_new") - self.load("m_old")
        y_k = self.load("g_new") - self.load("g_old")

        # Determine the shape of the memory map (length of model, length of mem)
        m = len(s_k)
        n = self.memory

        # Initial iteration, need to create the memory map
        if self.memory_used == 0:
            s = np.memmap(filename="LBFGS/S", mode="w+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename="LBFGS/Y", mode="w+", dtype="float32",
                          shape=(m, n))
            # Store the model and gradient differences in memmaps
            s[:, 0] = s_k
            y[:, 0] = y_k
            self.memory_used = 1
        # Subsequent iterations
        else:
            s = np.memmap(filename="LBFGS/S", mode="r+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename="LBFGS/Y", mode="r+", dtype="float32",
                          shape=(m, n))
            # Shift all stored memory by one index to make room for latest mem
            s[:, 1:] = s[:, :-1]
            y[:, 1:] = y[:, :-1]
            # Store the latest model and gradient in first index
            s[:, 0] = s_k
            y[:, 0] = y_k

            # Keep track of the memory used
            if self.memory_used < self.memory:
                self.memory_used += 1

        return s, y

    def apply(self, q, s=None, y=None):
        """
        Applies L-BFGS inverse Hessian to given vector

        :type q: np.array
        :param q: gradient direction too apply L-BFGS to
        :type s: np.memmap
        :param s: memory of model differences
        :type y: np.memmap
        :param y: memory of gradient direction differences
        :rtype r: np.array
        :return r: new search direction from application of L-BFGS
        """
        unix.cd(self.path)

        # If no memmaps are given as arguments, instantiate them
        if s is None or y is None:
            m = len(q)
            n = self.memory
            s = np.memmap(filename="LBFGS/S", mode="w+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename="LBFGS/Y", mode="w+", dtype="float32",
                          shape=(m, n))

        # First matrix product
        # Recursion step 2 from appendix A of Modrak & Tromp 2016
        kk = self.memory_used
        rh = np.zeros(kk)
        al = np.zeros(kk)
        for ii in range(kk):
            rh[ii] = 1 / np.dot(y[:, ii], s[:, ii])
            al[ii] = rh[ii] * np.dot(s[:,ii], q)
            q = q - al[ii] * y[:,ii]

        # Apply a preconditioner
        if self.precond:
            r = self.precond(q)
        else:
            r = q

        # Use scaling M3 proposed by Liu and Nocedal 1989
        sty = np.dot(y[:,0], s[:,0])
        yty = np.dot(y[:,0], y[:,0])
        r *= sty/yty

        # Second matrix product
        # Recursion step 4 from appendix A of Modrak & Tromp 2016
        for ii in range(kk - 1, -1, -1):
            be = rh[ii] * np.dot(y[:, ii], r)
            r = r + s[:, ii] * (al[ii] - be)

        return r

    def restart(self):
        """
        Discards history from memmaps and resets counters
        """
        self.iter = 1
        self.memory_used = 0

        unix.cd(self.path)
        s = np.memmap(filename="LBFGS/S", mode="r+")
        y = np.memmap(filename="LBFGS/Y", mode="r+")
        s[:] = 0.
        y[:] = 0.

    def check_status(self, g, r):
        """
        Check the status of the apply(), determine if restart necessary

        :type g: np.array
        :param g: current gradient direction
        :type r: np.array
        :param r: new gradient direction
        :rtype: int
        :return: status based on status check
        """
        theta = 180. * np.pi ** -1 * angle(g, r)

        if not 0. < theta < 90.:
            if self.verbose:
                print("restarting LBFGS... [not a descent direction]")
            return 1
        elif theta > 90. - self.thresh:
            if self.verbose:
                print("restarting LBFGS... [practical safeguard]")
            return 1
        else:
            return 0


