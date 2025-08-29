#!/usr/bin/env python3
"""
Nonlinear conjugate gradient method for optimization
"""
import numpy as np

from seisflows import logger
from seisflows.optimize.gradient import Gradient
from seisflows.tools.math import dot
from seisflows.tools.specfem_model import Model
from seisflows.tools import msg
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

        # Be honest!
        logger.warning(msg.mjr("THIS OPTIMIZATION MODULE IS NOT WELL TESTED "
                               "USE AT YOUR OWN RISK"))

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
        _acceptable_calc_beta = ["fletcher_reeves"]  # "pollak_ribere" 
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

        np.savez(file=self.path._checkpoint, **checkpoint_dict)  # NOQA

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
        g_new = Model(self.path._g_new)
        g_old = Model(self.path._g_old)
        p_old = Model(self.path._p_old)

        # CASE 1: If first iteration, search direction is the current gradient
        if self._NLCG_iter == 1:
            logger.info("first NLCG iteration, setting search direction "
                        "as inverse gradient")
            # p_new = -g
            return g_new.apply(actions=["*"], values=[-1], 
                               export_to=self.path._p_new)
        # CASE 2: Force restart if the iterations have surpassed the maximum
        # number of allowable iter
        elif self._NLCG_iter > self.NLCG_max:
            logger.info("restarting NLCG due to periodic restart condition. "
                        "setting search direction as inverse gradient")
            self.restart()
            # p_new = -g
            return g_new.apply(actions=["*"], values=[-1], 
                               export_to=self.path._p_new)
        # Normal NLCG direction computation
        else:
            if check_conjugacy(g_new, g_old) > self.NLCG_thresh:
                logger.info("restarting NLCG due to loss of conjugacy")
                self.restart()
                # p_new = -g
                return g_new.apply(actions=["*"], values=[-1], 
                                   export_to=self.path._p_new)

            # Compute beta, the scaling factor for the search direction. User
            # has chosen a specific underlying function to do this
            beta = self._calc_beta()

            # Optional step, apply preconditioner in place
            if self.preconditioner:
                g_new = self.precondition(q=g_new)

            # p_new = -1 * g + B * p_old (couldn't figure out how to do this in
            # one action so we apply three)
            p_new = g_new.apply(actions=["*"], values=[-1])
            p_old = p_old.apply(actions=["*"], values=[beta], 
                                export_to=self.path._scratch)
            p_new = p_new.apply(actions=["+"], values=[p_old])

            # Check restart conditions, return search direction and statusa
            if check_descent(p_new.get("vector"), g_new.get("vector")) > 0.:
                logger.info("restarting NLCG, not a descent direction")
                self.restart()
                # p_new = -g
                return g_new.apply(actions=["*"], values=[-1], 
                                   export_to=self.path._p_new)

        return p_new

    def restart(self):
        """
        Overwrite the Base restart class and include a restart of the NLCG
        """
        logger.info("restarting NLCG optimization algorithm")
        self._line_search.clear_search_history()
        self._restarted = True
        self._NLCG_iter = 1

    def _fletcher_reeves(self):
        """
        Calculating scale factor beta in the NLCG Algorithm following:
        Fletcher & Reeves, 1964

        :rtype: float
        :return: beta, the scale factor to apply to the old search direction to
            determine the new search direction
        """
        # Use current and old gradients
        g_new = Model(self.path._g_new)
        g_old = Model(self.path._g_old)

        if self.preconditioner:
            g_new_prcnd = self.precondition(g_new, export_to=self.path._scratch)
        else:
            g_new_prcnd = g_new
            
        return g_new_prcnd.dot(g_new) / g_old.dot(g_old)

    def _pollak_ribere(self):
        """
        Calculating scale factor beta in the NLCG Algorithm following:
        Polak & Ribiere, 1969

        :rtype: float
        :return: beta, the scale factor to apply to the old search direction to
            determine the new search direction
        """
        # BC Did not have time to finish updating this, need to figure out how
        # to store two scratch variables at the same time. See comment below
        raise NotImplementedError("This function is currently not working, if "
                                  "you would like to use it, please open a " 
                                  "GitHub issue")
    
        # Use current and old gradients
        g_new = Model(self.path._g_new)
        g_old = Model(self.path._g_old)
        p_old = Model(self.path._p_old)

        if self.preconditioner:
            g_new_prcnd = self.precondition(g_new, export_to=self.path._scratch)
        else:
            g_new_prcnd = g_new

        # This wont work because both `dg` and `g_new_prcnd` need to be stored
        # in scratch. I guess I could make another scratch directory or store it
        # somewhere else but I feel like there is a cleaner way
        dg = g_new.apply(action=["-"], values=g_old, 
                         export_to=self.path._scratch)
        
        return g_new_prcndt.dot(dg) / g_old.dot(g_old)
        

def check_conjugacy(g_new, g_old):
    """
    Check for conjugacy between two vectors

    :type g_new: seisflows.tools.specfem_model.Model
    :param g_new: new search direction
    :type g_old: seisflows.tools.specfem_model.Model
    :param g_old: old search direction
    :rtype: float
    :return: an element that proves conjugacy
    """
    return abs(g_new.dot(g_old) / g_new.dot(g_new))


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
    return p_new.dot(g_new) / g_new.dot(g_new)



