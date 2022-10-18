#!/usr/bin/env python3
"""
Nonlinear conjugate gradient method for optimization
"""
import numpy as np

from seisflows import logger
from seisflows.optimize.gradient import Gradient
from seisflows.tools.math import dot
from seisflows.plugins import line_search as line_search_dir


class NLCG(Gradient):
    """
    NLCG Optimization
    -----------------
    Nonlinear conjugate gradient method

    Parameters
    ----------
    :type nlcg_max: int
    :param nlcg_max: NLCG periodic restart interval, should be between 1
        and infinity
    :type nlcg_thresh: NLCG conjugacy restart threshold, should be
        between 1 and infinity
    :type calc_beta: str
    :param calc_beta: method to calculate the parameter 'beta' in the
        NLCG algorithm. Available: 'pollak_ribere', 'fletcher_reeves'

    Paths
    -----
    ***
    """
    __doc__ = Gradient.__doc__ + __doc__

    def __init__(self, nlcg_max=np.inf, nlcg_thresh=np.inf,
                 calc_beta="pollak_ribere", **kwargs):
        """NLCG-specific input parameters"""
        super().__init__(**kwargs)

        # Overwrite user-chosen line search. L-BFGS requires 'Backtrack'ing LS
        if self.line_search_method.title() != "Bracket":
            logger.warning(f"NLCG optimization requires 'bracket'ing line "
                           f"search. Overwritng '{self.line_search_method}'")
            self.line_search_method = "Bracket"
            self._line_search = getattr(
                line_search_dir, self.line_search_method)(
                step_count_max=self.step_count_max,
                step_len_max=self.step_len_max
            )

        self.NLCG_max = nlcg_max
        self.NLCG_thresh = nlcg_thresh

        # Check paramter validity
        _acceptable_calc_beta = ["pollak_ribere", "fletcher_reeves"]
        assert(calc_beta in _acceptable_calc_beta), (
            f"unacceptable `calc_beta`, must be in {_acceptable_calc_beta}"
        )
        self.calc_beta = calc_beta

        # Internally used parameters
        self._NLCG_iter = 0
        self._calc_beta = getattr(self, f"_{calc_beta}")

    def checkpoint(self):
        """
        Overwrite default checkpointing to store internal L-BFGS Attributes
        """
        super().checkpoint()

        checkpoint_dict = np.load(self.path._checkpoint)
        checkpoint_dict["NLCG_iter"] = self._NLCG_iter

        np.savez(file=self.path._checkpoint, **dict_out)  # NOQA

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

        :rtype: seisflows.tools.specfem.Model
        :return: search direction as a Model instance
        """
        self._NLCG_iter += 1

        # Load the current gradient direction
        g_new = self.load_vector("g_new")
        p_new = g_new.copy()

        # CASE 1: If first iteration, search direction is the current gradient
        if self._NLCG_iter == 1:
            logger.info("first NLCG iteration, setting search direction "
                        "as inverse gradient")
            p_new.update(vector=-1 * g_new.vector)
            restarted = False
        # CASE 2: Force restart if the iterations have surpassed the maximum
        # number of allowable iter
        elif self._NLCG_iter > self.NLCG_max:
            logger.info("restarting NLCG due to periodic restart "
                        "condition. setting search direction as inverse "
                        "gradient")
            self.restart()
            p_new.update(vector=-1 * g_new.vector)
            restarted = True
        # Normal NLCG direction compuitation
        else:
            # Compute search direction
            g_old = self.load_vector("g_old")
            p_old = self.load_vector("p_old")
            beta = self._calc_beta(g_new.vector, g_old.vector)

            # Apply preconditioner and calc. scale factor for search dir. (beta)
            _p_new_vec = (
                    -1 * self._precondition(g_new.vector) + beta * p_old.vector
            )
            p_new.update(vector=_p_new_vec)

            # Check restart conditions, return search direction and statusa
            if check_conjugacy(g_new.vector, g_old.vector) > self.NLCG_thresh:
                logger.info("restarting NLCG due to loss of conjugacy")
                self.restart()
                p_new.update(vector=-1 * g_new.vector)
                restarted = True
            elif check_descent(p_new.vector, g_new.vector) > 0.:
                logger.info("restarting NLCG, not a descent direction")
                self.restart()
                p_new.update(vector=-1 * g_new.vector)
                restarted = True
            else:
                # p_new = p_new
                restarted = False

        # Save values to disk and memory
        self._restarted = restarted

        return p_new

    def restart(self):
        """
        Overwrite the Base restart class and include a restart of the NLCG
        """
        logger.info("restarting NLCG optimization algorithm")

        g_new = self.load_vector("g_new")
        p_new = g_new.copy()
        p_new.update(vector=-1 * g_new.vector)
        self.save_vector("p_new", p_new)

        self._line_search.clear_search_history()
        self._restarted = 1
        self._NLCG_iter = 1

    def _fletcher_reeves(self, g_new, g_old):
        """
        One method for calculating beta in the NLCG Algorithm from
        Fletcher & Reeves, 1964

        :type g_new: np.array
        :param g_new: new search direction
        :type g_old: np.array
        :param g_old: old search direction
        :rtype: float
        :return: beta, the scale factor to apply to the old search direction to
            determine the new search direction
        """
        num = dot(self._precondition(g_new), g_new)
        den = dot(g_old, g_old)
        beta = num / den
        return beta

    def _pollak_ribere(self, g_new, g_old):
        """
        One method for calculating beta in the NLCG Algorithm from
        Polak & Ribiere, 1969

        :type g_new: np.array
        :param g_new: new search direction
        :type g_old: np.array
        :param g_old: old search direction
        :rtype: float
        :return: beta, the scale factor to apply to the old search direction to
            determine the new search direction
        """
        num = dot(self._precondition(g_new), g_new - g_old)
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



