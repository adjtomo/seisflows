#!/usr/bin/env python
"""
This is the custom class for an LBFGS optimization schema.
It supercedes the `seisflows.optimize.base` class
"""
import sys
import logging
import numpy as np

from seisflows3.config import custom_import, SeisFlowsPathsParameters
# from seisflows3.plugins import optimize

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class LBFGS(custom_import("optimize", "base")):
    """
    The Limited memory BFGS algorithm
    Calls upon seisflows.plugin.optimize.LBFGS to accomplish LBFGS algorithm
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type LBFGS: Class
        :param LBFGS: plugin LBFGS class that controls the machinery of the
            L-BFGS optimization schema
        :type LBFGS_iter: int
        :param LBFGS_iter: an internally used iteration that differs from
            optimization iter. Keeps track of internal LBFGS memory of previous
            gradients. If LBFGS is restarted, the LBFGS_iter iteration is reset,
            but the optization iteration.
        :type memory_used: int
        :param memory_used: bookkeeping to see how many previous
            gradients have been stored to internal memory. Should not exceed
            PAR.LBFGSMEM
        :type LBFGS_dir: str
        :param LBFGS_dir: location to store LBFGS internal memory
        :type y_file: str
        :param y_file: path to store memory of the gradient differences
            i.e., `g_new - g_old`
        :type s_file: str
        :param s_file: path to store memory of the model differences
            i.e., `m_new - m_old`
        """
        super().__init__()
        self.LBFGS_iter = 0
        self.memory_used = 0
        self.LBFGS_dir = "LBFGS"
        self.y_file = os.path.join(self.LBFGS_dir, "Y")
        self.s_file = os.path.join(self.LBFGS_dir, "S")


    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        # Define the Parameters required by this module
        sf.par("LINESEARCH", required=False, default="Backtrack", par_type=str,
               docstr="Algorithm to use for line search, see "
                      "seisflows.plugins.line_search for available choices")

        sf.par("LBFGSMEM", required=False, default=3, par_type=int,
               docstr="Max number of previous gradients to retain "
                      "in local memory")

        sf.par("LBFGSMAX", required=False, par_type=int, default="inf",
               docstr="LBFGS periodic restart interval, between 1 and 'inf'")

        sf.par("LBFGSTHRESH", required=False, default=0., par_type=float,
               docstr="LBFGS angle restart threshold")

        return sf

    def check(self, validate=True):
        """
        Checks parameters, paths, and dependencies
        """
        super().check(validate=False)
        if validate:
            self.required.validate()

        assert(PAR.LINESEARCH.upper() == "BACKTRACK"), \
            "LBFGS requires a Backtracking line search"

    def setup(self):
        """
        Set up the LBFGS optimization schema
        """
        super().setup()

        # Create a separate directory for LBFGS matters
        unix.cd(PATH.OPTIMIZE)
        unix.mkdir(self.LBFGS_dir)

    def compute_direction(self):
        """
        Call on the L-BFGS optimization machinery to compute a search
        direction using internally stored memory of previous gradients
        """
        self.logger.info(f"computing search direction with L-BFGS")
        self.LBFGS_iter += 1

        unix.cd(PATH.OPTIMIZE)

        # Load the current gradient direction, which is the L-BFGS search
        # direction if this is the first iteration
        g = self.load(self.g_new)

        # Restart condition or first iteration lead to setting search direction
        # as the inverse gradient (i.e., default to steepest descent)
        if self.LBFGS_iter == 1:
            self.logger.info("first L-BFGS iteration, setting search direction"
                             "as inverse gradient")
            self.save(self.p_new, -g)
            self.restarted = 0
            return
        elif self.LBFGS_iter > PAR.LBFGSMAX:
            self.logger.info("restarting L-BFGS due to periodic restart "
                             "condition. setting search direction as"
                             "inverse gradient")
            self.restart()
            self.save(self.p_new, -g)
            self.restarted = 1
            return

        # Update the search direction, apply the inverse Hessian
        # 'q' becomes the new search direction 'g'
        self.logger.info("applying inverse Hessian to gradient")
        s, y = self.update()
        q = self.apply(g, s, y)

        # Determine if the new search direction is appropriate by checking its
        # angle to the previous search direction
        status = self.check_status(g, q)

        if status == 0:
            self.logger.info("new L-BFGS search direction found")
            self.save(self.p_new, -q)
            self.restarted = status
        elif status == 1:
            self.logger.info("new search direction not appropriate, defaulting "
                             "to inverse gradient")
            self.restart()
            self.save(self.p_new, -g)
            self.restarted = status

    def restart(self):
        """
        Overwrite the Base restart class and include a restart of the LBFGS
        """
        super().restart()

        self.logger.info("restarting L-BFGS optimization algorithm by clearing "
                         "internal memory")
        self.LBFGS_iter = 1
        self.memory_used = 0

        unix.cd(self.path)
        s = np.memmap(filename=self.s_file, mode="r+")
        y = np.memmap(filename=self.y_file, mode="r+")
        s[:] = 0.
        y[:] = 0.

    def update(self):
        """
        Updates L-BFGS algorithm history

        .. note::
            Because models are large, and multiple iterations of models need to
            be stored in memory, previous models are stored as `memmaps`,
            which allow for access of small segments of large files on disk,
            without reading the entire file. Memmaps are array like objects.

        .. note::
            Notation for s and y taken from Liu & Nocedal 1989
            iterate notation: sk = x_k+1 - x_k and yk = g_k+1 - gk

        :rtype s: np.memmap
        :return s: memory of the model differences `m_new - m_old`
        :rtype y: np.memmap
        :return y: memory of the gradient differences `g_new - g_old`
        """
        unix.cd(PATH.OPTIMIZE)

        # Determine the iterates for model m and gradient g
        s_k = self.load(self.m_new) - self.load(self.m_old)
        y_k = self.load(self.g_new) - self.load(self.g_old)

        # Determine the shape of the memory map (length of model, length of mem)
        m = len(s_k)
        n = PAR.LBFGSMEM

        # Initial iteration, need to create the memory map
        if self.memory_used == 0:
            s = np.memmap(filename=self.s_file, mode="w+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename=self.y_file, mode="w+", dtype="float32",
                          shape=(m, n))
            # Store the model and gradient differences in memmaps
            s[:, 0] = s_k
            y[:, 0] = y_k
            self.memory_used = 1
        # Subsequent iterations will append to memory maps
        else:
            s = np.memmap(filename=self.s_file, mode="r+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename=self.y_file, mode="r+", dtype="float32",
                          shape=(m, n))
            # Shift all stored memory by one index to make room for latest mem
            s[:, 1:] = s[:, :-1]
            y[:, 1:] = y[:, :-1]
            # Store the latest model and gradient in first index
            s[:, 0] = s_k
            y[:, 0] = y_k

            # Keep track of the memory used
            if self.memory_used < PAR.LBFGSMEM:
                self.memory_used += 1

        return s, y

    def apply(self, q, s=None, y=None):
        """
        Applies L-BFGS inverse Hessian to given vector

        :type q: np.array
        :param q: gradient direction to apply L-BFGS to
        :type s: np.memmap
        :param s: memory of model differences
        :param s: memory of model differences
        :type y: np.memmap
        :param y: memory of gradient direction differences
        :rtype r: np.array
        :return r: new search direction from application of L-BFGS
        """
        unix.cd(PATH.OPTIMIZE)

        # If no memmaps are given as arguments, instantiate them
        if s is None or y is None:
            m = len(q)
            n = PAR.LBFGSMEM
            s = np.memmap(filename=self.s_file, mode="w+", dtype="float32",
                          shape=(m, n))
            y = np.memmap(filename=self.y_file, mode="w+", dtype="float32",
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

        # Apply a preconditioner if available
        if self.precond:
            r = self.precond(q)
        else:
            r = q

        # Use scaling M3 proposed by Liu and Nocedal 1989
        sty = np.dot(y[:, 0], s[:, 0])
        yty = np.dot(y[:, 0], y[:, 0])
        r *= sty/yty

        # Second matrix product
        # Recursion step 4 from appendix A of Modrak & Tromp 2016
        for ii in range(kk - 1, -1, -1):
            be = rh[ii] * np.dot(y[:, ii], r)
            r = r + s[:, ii] * (al[ii] - be)

        return r

    def check_status(self, g, r):
        """
        Check the status of the apply() function, determine if restart necessary
        Return of 1 means restart, return of 0 means good to go.

        :type g: np.array
        :param g: current gradient direction
        :type r: np.array
        :param r: new gradient direction
        :rtype: int
        :return: status based on status check. 0=good, 1=bad
        """
        theta = 180. * np.pi ** -1 * angle(g, r)
        self.logger.info(f"new search direction is {theta:.2f}deg from current")

        if not 0. < theta < 90.:
            self.logger.info("restarting L-BFGS, theta not a descent direction")
            return 1
        elif theta > 90. - PAR.LBFGSTHRESH:
            self.logger.info("restarting L-BFGS due to practical safeguard")
            return 1
        else:
            return 0