#!/usr/bin/env python
"""
This is the base class for seisflows.optimize
This class provides the core utilities for the Seisflows optimization schema.
"""
import os
import sys
import numpy as np

from seisflows.plugins import line_search, preconds
from seisflows.tools import msg, unix
from seisflows.tools.err import ParameterError
from seisflows.tools.tools import loadnpy, savenpy
from seisflows.tools.math import angle, poissons_ratio
from seisflows.tools.seismic import Writer

# seisflows.config objects 
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class Base(object):
    """
    Nonlinear optimization abstract base class

    Base class on top of which steepest descent, nonlinear conjugate, quasi-
    Newton and Newton methods can be implemented.  Includes methods for
    both search direction and line search.

    To reduce memory overhead, vectors are read from disk rather than passed
    from calling routines. For example, at the beginning of compute_direction
    the current gradient is  read from  'g_new' and the resulting search
    direction is written to 'p_new'. As the inversion progresses, other
    information is stored as well.

    Variables:
        m_new - current model
        m_old - previous model
        m_try - line search model
        f_new - current objective function value
        f_old - previous objective function value
        f_try - line search function value
        g_new - current gradient direction
        g_old - previous gradient direction
        p_new - current search direction
        p_old - previous search direction
    """
    def check(self):
        """
        Checks parameters, paths, and dependencies

        Note:
        The default numerical parameters defined below should work well for a
        range of applications without manual tuning. If the nonlinear
        optimization procedure stagnates, it may be due to issues involving
        data quality or the choice of data misfit, data processing, or
        regularization parameters.  Problems in any of these areas usually
        manifest themselves through stagnation of the nonlinear optimization
        algorithm.
        """
        # line search algorithm
        if "LINESEARCH" not in PAR:
            setattr(PAR, "LINESEARCH", "Bracket")
        
        # Preconditioner
        if "PRECOND" not in PAR:
            setattr(PAR, "PRECOND", None)

        # Maximum number of trial steps
        if "STEPCOUNTMAX" not in PAR:
            setattr(PAR, "STEPCOUNTMAX", 10)

        # Initial step length as fraction of current model
        if "STEPLENINIT" not in PAR:
            setattr(PAR, "STEPLENINIT", 0.05)

        # Maximum step length as fraction of current model
        if "STEPLENMAX" not in PAR:
            setattr(PAR, "STEPLENMAX", 0.5)

        # Location of temporary files
        if "OPTIMIZE" not in PATH:
            setattr(PATH, "OPTIMIZE", f"{PATH.SCRATCH}/optimize")

        # Assertions
        if "WORKDIR" not in PATH:
            raise ParameterError

        if PAR.OPTIMIZE in ['base']:
            print(msg.CompatibilityError1)
            sys.exit(-1)

        if PAR.LINESEARCH:
            assert PAR.LINESEARCH in dir(line_search)

        if PAR.PRECOND:
            assert PAR.PRECOND in dir(preconds)

        if PAR.STEPLENINIT:
            assert 0. < PAR.STEPLENINIT

        if PAR.STEPLENMAX:
            assert 0. < PAR.STEPLENMAX

        if PAR.STEPLENINIT and PAR.STEPLENMAX:
            assert PAR.STEPLENINIT < PAR.STEPLENMAX

    def setup(self):
        """
        Sets up nonlinear optimization machinery
        """
        # Prepare line search machinery
        self.line_search = getattr(line_search, PAR.LINESEARCH)(
            step_count_max=PAR.STEPCOUNTMAX,
            path=os.path.join(PATH.WORKDIR, "output.optim"),
            verbose=PAR.VERBOSE
        )

        # Prepare preconditioner
        if PAR.PRECOND:
            self.precond = getattr(preconds, PAR.PRECOND)()
        else:
            self.precond = None

        # Prepare output logs
        self.writer = Writer(path=os.path.join(PATH.WORKDIR, "output.stats"))

        # Prepare scratch directory and save initial model
        unix.mkdir(PATH.OPTIMIZE)
        if "MODEL_INIT" in PATH:
            solver = sys.modules['seisflows_solver']
            m_new = solver.merge(solver.load(PATH.MODEL_INIT))
            self.save("m_new", m_new)
            self.check_model_parameters(m_new, "m_new")

    def compute_direction(self):
        """
        Computes search direction

        Note:
            This function implements steepest descent, for other algorithms,
            simply overload this method
        """
        g_new = self.load('g_new')
        if self.precond:
            p_new = -self.precond(g_new)
        else:
            p_new = -g_new
        self.save('p_new', p_new)
    
    def check_model_parameters(self, m, tag):
        """
        Check to ensure that the model parameters fall within the guidelines 
        of the solver. Print off min/max model parameters for the User.
        
        :type m: np.array
        :param m: model to check parameters of 
        :type tag: str
        :param tag: tag of the model to be used for more specific error msgs
        """
        vp, vs, rho, poissons = None, None, None, None
        # Distribute the model parameters based on the type of inversion
        if PAR.MATERIALS == "Elastic":
            if PAR.DENSITY == "Constant":
                vp, vs = np.split(m, 2)
            else:
                # include density in the model vector
                vp, vs, rho = np.split(m, 3)
        else:
            vp = m

        # Check the Poisson's ratio based on Specfem3D upper/lower bounds
        if vp is not None and vs is not None:
            poissons = poissons_ratio(vp=vp, vs=vs)
            if (poissons.min() < -1) or (poissons.max() > 0.5):
                msg = f"""
                
                ERROR CHECKING MODEL PARAMETERS

                    The Poisson's ratio of '{tag}' is out of bounds with respect 
                    to the range defined by Specfem3D (-1, 0.5). This will cause 
                    an error in the process xgenerate_databases. The model 
                    bounds were found to be:

                    {poissons.min():.2f} < PR < {poissons.max():.2f}

                """
                print(msg)
                sys.exit(-1)

        # Tell the User min and max values of the updated model
        print("\t\tModel Parameters")
        msg = "\t\t\t{minval:.2f} <= {val} <= {maxval:.2f}"
        for val, tag in zip([vp, vs, rho, poissons], 
                            ["Vp", "Vs", "rho", "Poissons"]):
            if val is not None:
                print(msg.format(minval=val.min(), val=tag, maxval=val.max()))

    def initialize_search(self):
        """
        Determines first step length in line search
        """
        # Load in and calucate the necessary variables
        m = self.load('m_new')
        g = self.load('g_new')
        p = self.load('p_new')
        f = self.loadtxt('f_new')
        norm_m = max(abs(m))
        norm_p = max(abs(p))
        gtg = self.dot(g, g)
        gtp = self.dot(g, p)

        # Restart line search if necessary
        if self.restarted:
            self.line_search.clear_history()

        # Optional step length safeguard
        if PAR.STEPLENMAX:
            self.line_search.step_len_max = PAR.STEPLENMAX * norm_m / norm_p
            print(f"\t\tStep length max {self.line_search.step_len_max}")

        # Determine initial step length
        alpha, _ = self.line_search.initialize(0., f, gtg, gtp)

        # Optional initial step length override
        if PAR.STEPLENINIT and len(self.line_search.step_lens) <= 1:
            alpha = PAR.STEPLENINIT * norm_m / norm_p
            print("\t\tStep length override due to PAR.STEPLENINIT")

        # The new model is the old model, scaled by the step direction and
        # gradient threshold to remove any outlier values
        m_try = m + alpha * p

         # Write model corresponding to chosen step length
        self.save("m_try", m_try)
        self.savetxt("alpha", alpha)

        # Check the new model and update the User on a few parameters
        self.check_model_parameters(m_try, "m_try")

    def update_search(self):
        """
        Updates line search status and step length

        Status codes:
            status > 0  : finished
            status == 0 : not finished
            status < 0  : failed
        """
        alpha, status = self.line_search.update(self.loadtxt("alpha"),
                                                self.loadtxt("f_try"))

        # write model corresponding to chosen step length
        if status >= 0:
            m = self.load("m_new")
            p = self.load("p_new")
            self.savetxt("alpha", alpha)
            m_try = m + alpha * p
            self.save("m_try", m_try)
            self.check_model_parameters(m_try, "m_try")

        return status

    def finalize_search(self):
        """
        Prepares algorithm machinery and scratch directory for next model update

        Removes old model/search parameters, moves current parameters to old,
        sets up new current parameters and writes statistic outputs
        """
        # m = self.load('m_new')  # unusued variable
        g = self.load('g_new')
        p = self.load('p_new')
        x = self.line_search.search_history()[0]
        f = self.line_search.search_history()[1]

        # Clean scratch directory
        unix.cd(PATH.OPTIMIZE)
        # Remove the old model parameters
        if self.iter > 1:
            for fid in ["m_old", "f_old", "g_old", "p_old", "s_old"]:
                unix.rm(fid)

        # Rename current model parameters to "_old" for new search
        unix.mv("m_new", "m_old")
        unix.mv("f_new", "f_old")
        unix.mv("g_new", "g_old")
        unix.mv("p_new", "p_old")

        # Setup the current model parameters
        unix.mv("m_try", "m_new")
        self.savetxt("f_new", f.min())

        # Output latest statistics
        self.writer("factor",
                    -self.dot(g, g) ** -0.5 * (f[1] - f[0]) / (x[1] - x[0]))
        self.writer("gradient_norm_L1", np.linalg.norm(g, 1))
        self.writer("gradient_norm_L2", np.linalg.norm(g, 2))
        self.writer("misfit", f[0])
        self.writer("restarted", self.restarted)
        self.writer("slope", (f[1] - f[0]) / (x[1] - x[0]))
        self.writer("step_count", self.line_search.step_count)
        self.writer("step_length", x[f.argmin()])
        self.writer("theta", 180. * np.pi ** -1 * angle(p, -g))

        self.line_search.writer.newline()

    def retry_status(self):
        """
        After a failed line search, this determines if restart is worthwhile
        by checking, in effect, if the search direction was the same as gradient
        direction
        """
        g = self.load("g_new")
        p = self.load("p_new")
        theta = angle(p,-g)

        if PAR.VERBOSE:
            print(f" theta: {theta:6.3f}")

        thresh = 1.e-3
        if abs(theta) < thresh:
            return 0
        else:
            return 1

    def restart(self):
        """
        Restarts nonlinear optimization algorithm.
        Keeps current position in model space, but discards history of
        nonlinear optimization algorithm in an attempt to recover from
        numerical stagnation.
        """
        g = self.load("g_new")
        self.save("p_new", -g)
        self.line_search.clear_history()
        self.restarted = 1
        self.line_search.writer.iter -= 1
        self.line_search.writer.newline()

    def dot(self, x, y):
        """
        Utility function to computes inner product between vectors

        :type x: np.array
        :param x: vector 1
        :type y: np.array
        :param y: vector 2
        """
        return np.dot(np.squeeze(x), np.squeeze(y))

    def load(self, filename):
        """
        Reads vectors from disk

        :type filename: str
        :param filename: filename to read from
        :return:
        """
        return loadnpy(os.path.join(PATH.OPTIMIZE, filename))

    def save(self, filename, array):
        """
        Writes vectors to disk

        :type filename: str
        :param filename: filename to read from
        :type array: np.array
        :param array: array to be saved
        :return:
        """
        savenpy(os.path.join(PATH.OPTIMIZE, filename), array)

    def loadtxt(self, filename):
        """
        Reads scalars from disk

        :type filename: str
        :param filename: filename to read from
        :type array: np.array
        :param array: array to be saved
        :return:
        """
        return float(np.loadtxt(os.path.join(PATH.OPTIMIZE, filename)))

    def savetxt(self, filename, scalar):
        """
        Writes scalars to disk

        :type filename: str
        :param filename: filename to read from
        :type scalar: float
        :param scalar: value to write to disk
        :return:
        """
        np.savetxt(os.path.join(PATH.OPTIMIZE, filename), [scalar], "%11.6e")


