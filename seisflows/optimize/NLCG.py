#!/usr/bin/env python3
"""
This is the custom class for an NLCG optimization schema.
It inherits from the `seisflows.optimize.gradient.Gradient` class
"""
from seisflows.optimize.gradient import Gradient
from seisflows.tools import unix
from seisflows.tools.math import dot


class NLCG(Gradient):
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
    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        :type NLCG_iter: Class
        :param NLCG_iter: an internally used iteration that differs from
            optimization iter. Keeps track of internal NLCG memory.
        """
        super().__init__()

        self.required.par(
            "NLCGMAX", required=False, default="null", par_type=float,
            docstr="NLCG periodic restart interval, between 1 and inf"
        )
        self.required.par(
            "NLCGTHRESH", required=False, default="null", par_type=float,
            docstr="NLCG conjugacy restart threshold, between 1 and inf"
        )
        self.NLCG_iter = 0
        self.calc_beta = pollak_ribere  # !!! Allow the user to choose this fx?

    def check(self, validate=True):
        """
        Checks parameters, paths, and dependencies
        """
        super().check(validate=validate)

        assert(self.par.LINESEARCH.upper() == "BRACKET"), \
            f"NLCG requires a bracketing line search algorithm"

    def setup(self):
        """Inherit from optimize.gradient.Gradient"""
        self.setup()

    def finalize(self):
        """Inherit from optimize.gradient.Gradient"""
        self.finalize()

    def load(self, name):
        """Inherit from optimize.gradient.Gradient"""
        return self.load(name=name)

    def save(self, name, vector):
        """Inherit from optimize.gradient.Gradient"""
        self.save(name=name, vector=vector)

    def compute_direction(self):
        """
        Compute search direction using the Nonlinear Conjugate Gradient method
        The potential outcomes when computing direction with NLCG

        1. First iteration of an NLCG optimization, search direction is
            the inverse gradient
        2. NLCG internal iteration ticks over the maximum allowable number of
            iterations, force a restart condition, search direction is the
            inverse gradient
        3. New NLCG search direction does not have conjugacy with previous
            search direction, force restart, inverse gradient search direction
        4. New NLCG search direction is not a descent direction,
            force restart, inverse gradient search direction
        5. New NLCG search direction has conjugacy and is a descent direction
            and is set as the new search direction.
        """
        self.logger.debug(f"computing search direction with NLCG")
        self.NLCG_iter += 1

        unix.cd(self.path.OPTIMIZE)

        # Load the current gradient direction
        g_new = self.load("g_new")

        # CASE 1: If first iteration, search direction is the current gradient
        if self.NLCG_iter == 1:
            self.logger.info("first NLCG iteration, setting search direction"
                             "as inverse gradient")
            p_new = -g_new
            restarted = 0
        # CASE 2: Force restart if the iterations have surpassed the maximum
        # number of allowable iter
        elif self.NLCG_iter > self.par.NLCGMAX:
            self.logger.info("restarting NLCG due to periodic restart "
                             "condition. setting search direction as inverse "
                             "gradient")
            self.restart()
            p_new = -g_new
            restarted = 1
        # Normal NLCG direction compuitation
        else:
            # Compute search direction
            g_old = self.load("g_old")
            p_old = self.load("p_old")

            # Apply preconditioner and calc. scale factor for search dir. (beta)
            if self.precond:
                beta = self.calc_beta(g_new, g_old, self.precond)
                p_new = -self.precond(g_new) + beta * p_old
            else:
                beta = self.calc_beta(g_new, g_old)
                p_new = -g_new + beta * p_old

            # Check restart conditions, return search direction and status
            if check_conjugacy(g_new, g_old) > self.par.NLCGTHRESH:
                self.logger.info("restarting NLCG due to loss of conjugacy")
                self.restart()
                p_new = -g_new
                restarted = 1
            elif check_descent(p_new, g_new) > 0.:
                self.logger.info("restarting NLCG, not a descent direction")
                self.restart()
                p_new = -g_new
                restarted = 1
            else:
                p_new = p_new
                restarted = 0

        # Save values to disk and memory
        self.save("p_new", p_new)
        self.restarted = restarted

    def restart(self):
        """
        Overwrite the Base restart class and include a restart of the NLCG
        """
        super().restart()
        self.NLCG_iter = 1

    def write_stats(self):
        """Inherit from optimize.gradient.Gradient"""
        self.write_stats()

    def check_model(self, m, min_pr=-1, max_pr=0.5):
        """Inherit from optimize.gradient.Gradient"""
        self.check_model(m=m, min_pr=min_pr, max_pr=max_pr)


def fletcher_reeves(g_new, g_old, precond=lambda x: x):
    """
    One method for calculating beta in the NLCG Algorithm from
    Fletcher & Reeves, 1964

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
    One method for calculating beta in the NLCG Algorithm from
    Polak & Ribiere, 1969

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
    num = dot(precond(g_new), g_new - g_old)
    den = dot(g_old, g_old)
    beta = num / den
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



