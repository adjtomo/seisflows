#!/usr/bin/env python3
"""
L-BFGS (Limited memory Broyden–Fletcher–Goldfarb–Shanno) algorithm for solving
nonlinear optimization problems.

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

TODO store LBFGS_iter during checkpointing
"""
import os
import numpy as np

from seisflows import logger
from seisflows.optimize.gradient import Gradient
from seisflows.tools import unix
from seisflows.tools.msg import DEG
from seisflows.tools.math import angle
from seisflows.tools.specfem_model import Model
from seisflows.plugins import line_search as line_search_dir


class LBFGS(Gradient):
    """
    L-BFGS Optimization
    -------------------
    Limited memory BFGS nonlinear optimization algorithm

    Parameters
    ----------
    :type lbfgs_mem: int
    :param lbfgs_mem: L-BFGS memory. Max number of previous gradients to
        retain in local memory for approximating the objective function.
    :type lbfgs_max: L-BFGS periodic restart interval. Must be
        1 <= lbfgs_max <= infinity.
    :type lbfgs_thresh: L-BFGS angle restart threshold. If the angle between
        the current and previous search direction exceeds this value,
        optimization algorithm will be restarted.

    Paths
    -----
    ***
    """
    __doc__ = Gradient.__doc__ + __doc__

    def __init__(self, lbfgs_mem=3, lbfgs_max=np.inf, lbfgs_thresh=0.,
                 **kwargs):
        """Instantiate L-BFGS specific parameters"""
        super().__init__(**kwargs)

        # Overwrite user-chosen line search. L-BFGS requires 'Backtrack'ing LS
        if self.line_search_method.title() != "Backtrack":
            self.line_search_method = "Backtrack"
            self._line_search = getattr(
                line_search_dir, self.line_search_method)(
                    step_count_max=self.step_count_max,
            )

        self.LBFGS_mem = lbfgs_mem
        self.LBFGS_max = lbfgs_max
        self.LBFGS_thresh = lbfgs_thresh

        # Set new L-BFGS dependent paths for storing previous gradients
        self.path._LBFGS = os.path.join(self.path.scratch, "LBFGS")

        # `r` defines the inverse Hessian used for search direction compute
        # in `compute_direction()`
        self.path["_r"] = os.path.join(self.path._LBFGS, "r")

        # These are used to store LBFGS memory which include model and gradient
        # differences
        for i in range(self.LBFGS_mem):
            # `y` = model difference (m_i+1 - m_i)
            self.path[f"_y{i}"] = os.path.join(self.path._LBFGS, f"y{i}")

            # `s` = gradient difference (g_i+1 - g_i)
            self.path[f"_s{i}"] = os.path.join(self.path._LBFGS, f"s{i}")

        # Internally used memory parameters for the L-BFGS optimization algo.
        self._LBFGS_iter = 0
        self._memory_used = 0

    def setup(self):
        """
        Set up the LBFGS optimization schema
        """
        super().setup()
        unix.mkdir(self.path._LBFGS)

    def checkpoint(self):
        """
        Overwrite default checkpointing to store internal L-BFGS Attributes
        """
        super().checkpoint()
        checkpoint_dict = dict(np.load(self.path._checkpoint))
        checkpoint_dict["LBFGS_iter"] = self._LBFGS_iter
        checkpoint_dict["memory_used"] = self._memory_used

        np.savez(file=self.path._checkpoint, **checkpoint_dict)  # NOQA

    def load_checkpoint(self):
        """
        Counterpart to `optimize.checkpoint`. Loads a checkpointed optimization
        module from disk and sets internal tracking attributes. Adds additional
        functionality to restore internal L-BFGS attributes
        """
        super().load_checkpoint()

        # NumPy appends '.npz' when saving. Make sure we honor that.
        if not self.path._checkpoint.endswith(".npz"):
            fid = f"{self.path._checkpoint}.npz"
        else:
            fid = self.path._checkpoint

        if os.path.exists(fid):
            dict_in = np.load(file=fid)
            self._LBFGS_iter = int(dict_in["LBFGS_iter"])
            self._memory_used = int(dict_in["memory_used"])

    def compute_direction(self):
        """
        Call on the L-BFGS optimization machinery to compute a search
        direction using internally stored memory of previous gradients.
        
        Exports new search direction to `p_new`

        The potential outcomes when computing direction with L-BFGS:
        1. First iteration of L-BFGS optimization, search direction is defined
            as the inverse gradient
        2. L-BFGS internal iteration ticks over the maximum allowable number of
            iterations, force a restart condition, search direction is the
            inverse gradient
        3. New search direction vector is too far from previous direction,
            force a restart, search direction is inverse gradient
        4. New search direction is acceptably angled from previous,
            becomes the new search direction
        """
        self._LBFGS_iter += 1
        logger.info(f"LBFGS iteration incremented -> {self._LBFGS_iter}")

        # Load the current gradient direction 
        g = Model(path=self.path._g_new)  

        # First LBFGS iteration means we default to gradient descent p = -g
        if self._LBFGS_iter == 1:
            logger.info("first L-BFGS iteration, default to gradient descent "
                        "(P = -G)")
            # Search direction `p_new` is negative gradient 
            p_new = g.apply(actions=["*"], values=[-1], 
                            export_to=self.path._p_new)
        # Optional periodic restart condition to recalibrate w/ gradient desc.
        elif self._LBFGS_iter > self.LBFGS_max:
            logger.info("restarting L-BFGS due to periodic restart condition. "
                        "setting search direction as inverse gradient (P = -G)")
            self.restart()

            # New search direction `p_new` is negative gradient
            p_new = g.apply(actions=["*"], values=[-1], 
                            export_to=self.path._p_new)
        # Normal LBFGS direction computation
        else:
            # Update the search direction, apply the inverse Hessian such that
            # 'q' becomes the new search direction 'g'
            logger.info("applying inverse Hessian to gradient")
            r = self._apply_inverse_hessian_to_gradient(q=g)  # new grad. `r`

            # Check the status of the new search direction by calculating the
            # angle between the new gradient `r` and the old gradient `g`
            theta = (180. / np.pi) * g.angle(r)  # units: deg
            logger.info(f"new search direction `r` {theta:.2f}{DEG} from `g`")

            # Determine if the new search direction is appropriate by checking
            # its angle to the previous search direction
            if 0. < theta < (90. - self.LBFGS_thresh):
                logger.info("new L-BFGS direction `r` is a descent direction")
                # New search direction `p_new` is p = -r
                p_new = r.apply(actions=["*"], values=[-1], 
                                export_to=self.path._p_new)
            else:
                logger.info("new search direction not a descent direction, " 
                            "falling back to gradient descent")
                self.restart()

                if theta > (90. - self.LBFGS_thres):
                    logger.debug(f"restarting L-BFGS due to `LBFGS_thresh`: "
                                 f"angle ({theta:.2f}) > 90 - "
                                 f"safeguard ({self.LBFGS_thresh})")
                    
                # New search direction `p_new` is negative gradient
                p_new = g.apply(actions=["*"], values=[-1], 
                                export_to=self.path._p_new)

        return p_new

    def restart(self):
        """
        Restart the L-BFGS optimization algorithm by clearing out stored
        gradient memory. The workflow will need to call `compute_direction` and
        `initialize_search` afterwards to properly re-instantiate the line
        search machinery. 
        """
        logger.info("restarting L-BFGS optimization algorithm")

        # Clear internal memory back to default values
        self._line_search.clear_search_history()
        self._restarted = True
        self._LBFGS_iter = 0  # will be incremented to 1 by `compute_direction`
        self._memory_used = 0

        # Clear out previous model difference (s_k) and 
        # gradient difference (y_k) information 
        for i in self.LBFGS_mem:
            Model(os.path.join(self.path[f"_y{i}"])).clear()
            Model(os.path.join(self.path[f"_s{i}"])).clear()

    def _apply_inverse_hessian_to_gradient(self, q):
        """
        Applies L-BFGS inverse Hessian to given vector following Appendix A
        Recursion from Modrak and Tromp (2016). Ultimately calculates a new
        gradient vector `r based on L-BFGS Algorithm. Does not return but 
        updates path `LBFGS._r` where the updated gradient vector is stored.

        In words the Recursion step is as follows:
        (1) Set q = gk
            i = k− 1
            j= min (k, l) 
        (2) Perform j times: 
            lambda_i = s_i^T * q / y_i^T * s_i, 
            q -= lambda_i * y_i
            i -= 1
        (3) Set gamma = s^T_(k-1) * y_(k-1) / y^T_(k-1) * y_k-1
            r = gamma * D * q
            i = k - j
        (4) Perform j times
            mu = y_i^T * r / y_i^T
            r += s_i(lambda_i - mu)
            i += 1
        (5) End with result: p_k = -r

        :type q: seisflows.tools.specfem_model.Model
        :param q: current gradient, stored in `g_new`
        """
        # Step 0: Update memory for L-BFGS
        self._update_search_history()

        # Step 1: Set q=g_k, i=k-1, j=min(k,l)
        k = self._memory_used
        j = min(k, self.LBFGS_mem)

        # Used to store dot products so we don't have to recompute
        s_dot_q = np.zeros(k)
        y_dot_s = np.zeros(k)
        
        # Step 2: Perform j timesz 
        # Note that we are storing intermediate variables in the `scratch` path
        # because they do not need to be retained
        for i in range(j):
            y_i = Model(self.path[f"_y{i}"])
            s_i = Model(self.path[f"_s{i}"])

            # lambda_i = s_i^T * q / y_i^T * s_i
            s_dot_q[i] = s_i.dot(q)    
            y_dot_s[i] = y_i.dot(s_i)  # former: rh = 1 / y_dot_s
            lambda_i = s_dot_q[i] / y_dot_s[i]  # former al = s_dot_q / y_dot_s

            # q -= lambda_i * y_i
            # Note the math is rearranged so we can make this slightly more
            # efficient but it's solving the same operation. Export to scratch
            q = y_i.apply(actions=["*", "*", "+"], values=[lambda_i, -1, q],
                          export_to=self.path._scratch)

        # Optional step, apply preconditioner in place on the scratch path
        if self.preconditioner:
            r = self.precondition(q=q)
        else:
            r = q
            
        # Step 3: Use scaling M3 proposed by Liu & Nocedal 1989
        y_0 = Model(self.path._y0)
        s_0 = Model(self.path._s0)
        sty_over_yty = s_0.dot(y_0) / y_0.dot(y_0)

        # We now use the actual `r` path that will provide the final quantity
        # `r` which defines the scaled gradient
        r = r.apply(actions=["*"], values=[sty_over_yty], 
                    export_to=self.path._r)
        
        # Clear out the scratch directory now that we don't need it anymore
        q.clear()

        # Step 4: Second matrix product, calculate scaled gradient `r`
        for i in range(k - 1, -1 ,-1):  
            # Load vectors from disk
            y_i = Model(self.path[f"_y{i}"])
            s_i = Model(self.path[f"_s{i}"])

            # mu = y_i^T * r / y_i^T * s_i
            mu = y_i.dot(r) / y_dot_s[i]  # previously defined as `be`
            lambda_i = s_dot_q[i] / y_dot_s[i]

            # r += s_i(lambda_i - mu)
            # Again we slightly rearrange the math to make the operation more
            # efficient but in the end we update the `r` model vector
            r = s_i.apply(actions=["*", "+"], values=[lambda_i - mu, r],
                          export_to=self.path._r)
            
        return Model(self.path._r)

    def _update_search_history(self):
        """
        Updates L-BFGS algorithm history by removing old memory and storing
        new model differences (s_k) and gradient differences (y_k)

        .. note::
            This algorithm is notated in Modrak & Tromp 2016 Appendix A
        
        .. note::
            Notation for s and y taken from Liu & Nocedal 1989
            iterate notation: sk = x_k+1 - x_k and yk = g_k+1 - gk
        """
        # If we have memory assigned, we need to make space for the latest entry
        # Essentially were shifting every index down by one (i -> i+1) and 
        # removing the last index to make way for new 0th index memory. In 
        # practice we're just renaming directories
        if self._memory_used > 0:
            if self._memory_used == self.LBFGS_mem:
                # Remove the oldest memory if we exceed the memory threshold
                unix.rm(self.path[f"_y{self._memory_used - 1}"])
                unix.rm(self.path[f"_s{self._memory_used - 1}"])

            # Iterate backwards, so for mem==3 we go: 2, 1 (ignore 0)
            for i in np.arange(self._memory_used, 0, -1):
                # Reassigning the next index to the current. E.g., if mem==3 we
                # start at i==2 which we don't need anymore, and move i==1 into 
                # the 2 spot.
                unix.mv(src=self.path[f"_y{i-1}"], dst=self.path[f"_y{i}"])
                unix.mv(src=self.path[f"_s{i-1}"], dst=self.path[f"_s{i}"])
        
        # Assign the current model and gradient differences to 0th index
        # s_k = m_new - m_old
        Model(self.path._m_new).apply(
            actions=["-"], values=[Model(self.path._m_old)],
            export_to=self.path._s0
            )
        # y_k = g_new - g_old
        Model(self.path._g_new).apply(
            actions=["-"], values=[Model(self.path._g_old)],
            export_to=self.path._y0
            )

        # Increment memory up to LBFGS memory threshold
        if self._memory_used < self.LBFGS_mem:
            self._memory_used += 1

    def _check_status(self, g, r):
        """
        Check the status of the apply() function, determine if restart necessary
        Return of False means restart, return of True means good to go.

        :type g: seisflows.tools.specfem_model.Model
        :param g: current gradient direction
        :type r: seisflows.tools.specfem_model.Model
        :param r: new gradient direction
        :rtype: bool
        :return: okay status based on status check (False==bad, True==good)
        """
        theta = 180. * np.pi ** -1 * g.angle(r)
        logger.info(f"new search direction: {theta:.2f}{DEG} from current")

        if not 0. < theta < 90.:
            logger.info("restarting L-BFGS, theta not a descent direction")
            return False
        elif theta > 90. - self.LBFGS_thresh:
            logger.info("restarting L-BFGS due to practical safeguard")
            return False
        else:
            return True

