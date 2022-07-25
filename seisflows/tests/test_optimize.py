"""
Test the optimization module by setting up a Rosenbrock optimization problem
and running line search
"""
import os
import pytest
import numpy as np
from seisflows.tools.specfem import Model
from seisflows.optimize.gradient import Gradient
from seisflows.optimize.LBFGS import LBFGS
from seisflows.optimize.NLCG import NLCG


@pytest.fixture
def rosenbrock():
    """
    Rosenbrock test problem for optimization library testing

    https://en.wikipedia.org/wiki/Rosenbrock_function
    """
    model_init = np.array([-1.2, 1])  # This is the guess for the global min
    model_true = np.array([1, 1])  # This is the actual minimum

    def objective_function(x):
        """
        Rosenbrock objective function which is defined mathematically as:

        f(x,y) = (a-x)^2 + b(y-x^2)^2

        where the global minimum is at (x,y) == (a, a^2)
        and typical constant values are: a==1, b==100
        """
        return np.array([((1 - x[0]) ** 2 + 100 * (-x[0] ** 2 + x[1]) ** 2)])

    def gradient(x):
        """
        Gradient of the objective function for Rosenbrock test
        """
        return np.array([-2*(1-x[0]) - 400*x[0]*(-x[0]**2+x[1]),
                         200*(- x[0]**2+x[1])])

    return model_init, model_true, objective_function, gradient


def rosenbrock_n(n=1E5):
    """
    N dimensional Rosenbrock test problem for optimization library testing

    https://en.wikipedia.org/wiki/Rosenbrock_function
    """
    model_init = 0.1 * np.ones(int(n))  # This is a guess for the global min
    model_true = np.ones(int(n))

    def objective_function(x):
        """
        Rosenbrock objective function
        """
        return sum(100 * (x[:-1]**2. - x[1:])**2. + (x[:-1] - 1.)**2)

    def gradient(x):
        """
        Gradient of the objective function for Rosenbrock test
        """
        g = np.zeros(int(n))
        g[1:-1] = -200 * (x[:-2] ** 2. - x[1:-1]) + \
                  400. * x[1:-1] * (x[1:-1] ** 2. - x[2:]) + \
                  2. * (x[1:-1] - 1.)

        g[0] = 400. * x[0] * (x[0] ** 2. - x[1]) + \
               2. * (x[0] - 1)

        g[-1] = -200. * (x[-2] ** 2. - x[-1])

        return g

    return model_init, model_true, objective_function, gradient


def test_gradient_descent_w_bracket(tmpdir, rosenbrock):
    """
    Test Gradient class with Rosenbrock problem
    """
    m_new, m_true, evaluate_objective_function, evaluate_gradient = rosenbrock()

    optimize = Gradient(workdir=tmpdir, line_search_method="bracket")
    optimize.setup()
    optimize.check()

    # Set up the optimization problem
    optimize.save_vector(name="m_new", m=m_new)
    optimize.initialize_search()

    while True:
        f_try = evaluate_objective_function(x=m_new)
        optimize.save_vector(name="f_try", m=f_try)

        g_new = evaluate_gradient(x=m_new)
        optimize.save_vector(name="g_new", m=g_new)

        optimize.
