"""
Test the optimization module by setting up a Rosenbrock optimization problem
and running line search
"""
import os
import pytest
import numpy as np
import matplotlib.pyplot as plt
from seisflows.tools.config import Dict
from seisflows.tools.model import Model
from seisflows.tools.math import angle
from seisflows.optimize.gradient import Gradient
from seisflows.optimize.LBFGS import LBFGS
from seisflows.optimize.NLCG import NLCG


def rosenbrock_objective_function(x):
    """
    Rosenbrock objective function used to test the optimization module.
    The Rosenbrock is defined mathematically as:

    f(x,y) = (a-x)^2 + b(y-x^2)^2

    where the global minimum is at (x,y) == (a, a^2)
    and typical constant values are: a==1, b==100

    https://en.wikipedia.org/wiki/Rosenbrock_function

    :type x: np.array
    :param x: input model [x, y] to feed into the objective function
    :rtype: float
    :return: misfit value for the given input model
    """
    return (1 - x[0]) ** 2 + 100 * (-x[0] ** 2 + x[1]) ** 2


def rosenbrock_gradient(x):
    """
    Define the gradient of the Rosenbrock function, i.e., the gradient of the
    `rosenbrock_objective_function`

    :type x: np.array
    :param x: input model [x, y] to feed into the gradient of the obj function
    :rtype: np.array
    :return: gradient vector of the Rosenbrock function
    """
    return np.array([-2 * (1 - x[0]) - 400 * x[0] * (-x[0] ** 2 + x[1]),
                     200 * (- x[0] ** 2 + x[1])])


@pytest.fixture
def setup_optimization_vectors(tmpdir):
    """
    Create vectors required for a Rosenbrock optimization problem to be used
    to test the line search and optimization algorithms.

    .. note::
        The optimization module requires a model (m_new), corresponding
        misfit of that model (f_new), the gradient of that misfit
        function (g_new) and a search direction (p_new; calculated by each
        individual optimization sub-class)
    """
    m_init = Model()
    m_init.model = Dict(x=[np.array([-1.2, 1])])  # Starting guess for Rosenbrock
    m_init.save(path=os.path.join(tmpdir, "m_new.npz"))

    # Calculate the misfit and save as a text file
    misfit = rosenbrock_objective_function(x=m_init.vector)
    np.savetxt(os.path.join(tmpdir, "f_new.txt"), [misfit])

    # Calcualte the gradient and save as a Model
    gradient = rosenbrock_gradient(x=m_init.vector)
    g_new = m_init.copy()
    g_new.update(vector=gradient)
    g_new.save(path=os.path.join(tmpdir, "g_new.npz"))


def test_gradient_compute_direction(tmpdir, setup_optimization_vectors):
    """
    Ensure that gradient computes the correct search direction. The gradient
    vector created by the Rosenbrock function is expected to be [-215.6, -88.]
    """
    optimize = Gradient(path_optimize=tmpdir)

    # Check that setup created the correct gradient vector
    g_new = optimize.load_vector("g_new")
    assert(g_new.vector.sum() == pytest.approx(-303.59, 3))

    # Gradient descent search direction is just the negative gradient
    p_new = optimize.compute_direction()
    assert(p_new.vector.sum() == pytest.approx(303.59, 3))


def test_gradient_initialize_search(tmpdir, setup_optimization_vectors):
    """
    Test that the optimization module properly initializes the line search by
    providing a new model and misfit value.

    Initialize search requres 4 vectors to properly intialize the search,
    they are 'm_new' (model), 'g_new' (gradient), 'p_new' (search direciton_
    and 'f_new' (misfit)
    """
    optimize = Gradient(path_optimize=tmpdir)
    p_new = optimize.compute_direction()
    p_new.save(path=os.path.join(tmpdir, "p_new"))

    m_try, alpha = optimize.initialize_search()

    assert(m_try.vector[0] == pytest.approx(-1.14, 3))
    assert(m_try.vector[1] == pytest.approx(1.02, 3))
    assert(alpha == pytest.approx(2.789e-4, 4))


def test_gradient_update_line_search(tmpdir, setup_optimization_vectors):
    """
    Ensure that updating the line search works as advertised, i.e., we get a
    status on how the line search should proceed
    """
    optimize = Gradient(path_optimize=tmpdir)
    p_new = optimize.compute_direction()
    p_new.save(path=os.path.join(tmpdir, "p_new"))

    m_try, alpha = optimize.initialize_search()
    optimize.save_vector("m_try", m_try)
    optimize.save_vector("alpha", alpha)

    f_try = rosenbrock_objective_function(m_try.vector)
    np.savetxt(os.path.join(tmpdir, "f_try.txt"), [f_try])

    m_try, status, status = optimize.update_line_search()
    assert(status == "TRY")
    assert(optimize._line_search.step_count == 1)
    assert(optimize._line_search.func_vals[1] == f_try)


def test_bracket_line_search(tmpdir, setup_optimization_vectors):
    """
    Run a small optimization problem to try to reduce the Rosenbrock
    objective function using the Bracket'ing line search method. Checks that the
    line search only takes a few steps and that the reduced misfit value is
    as expected.
    """
    optimize = Gradient(path_optimize=tmpdir, path_output=tmpdir,
                        line_search_method="bracket", step_count_max=100)

    # Make sure the initial misfit is high
    assert(optimize.load_vector("f_new") == pytest.approx(24.2, 1))

    # Calculate the initial search direction which is just the inverse gradient
    p_new = optimize.compute_direction()
    p_new.save(path=os.path.join(tmpdir, "p_new"))

    # Saves trial model 'm_try' and corresponding step length 'alpha'
    m_try, alpha = optimize.initialize_search()
    optimize.save_vector("m_try", m_try)
    optimize.save_vector("alpha", alpha)

    def line_search():
        """Each line search evaluation requires calculating a new model and
        corresponding misfit value"""
        m_try = optimize.load_vector("m_try")
        f_try = rosenbrock_objective_function(m_try.vector)
        np.savetxt(os.path.join(tmpdir, "f_try.txt"), [f_try])

        m_try, alpha, status = optimize.update_line_search()
        optimize.save_vector("m_try", m_try)
        optimize.save_vector("alpha", alpha)
        return status

    # Run a line search until acceptable misfit reduction
    while True:
        status = line_search()
        if status == "PASS":
            break

    assert(status == "PASS")  # pass
    assert(optimize.step_count == 4)  # Took 4 steps to reduce misfit
    # Make sure we have reduced the final misfit
    assert(min(optimize._line_search.func_vals) == pytest.approx(4.22, 3))

    # Change model names and reset line search
    optimize.finalize_search()

    # Check that the final model
    m_new = optimize.load_vector("m_new")
    m_old = optimize.load_vector("m_old")
    m_angle = angle(m_new.vector, m_old.vector)
    assert(m_angle == pytest.approx(0.3, 2))


def test_optimize_checkpoint_reload(tmpdir):
    """
    Checkpointing is used to store the status of the optimization module
    in the case of a failed or re-started workflow. Test that this works as
    advertised
    """
    rand_val = 123
    optimize = Gradient(path_optimize=tmpdir)
    optimize._line_search.step_count = rand_val
    optimize.checkpoint()

    new_optimize = Gradient(path_optimize=tmpdir)
    new_optimize.load_checkpoint()
    assert(new_optimize.step_count == rand_val)


def test_optimize_attempt_line_search_restart(tmpdir,
                                              setup_optimization_vectors):
    """
    Make sure we can tell when we're supposed to attempt a line search restart,
    i.e., when the gradient and search direction are the same
    """
    optimize = Gradient(path_optimize=tmpdir)
    p_new = optimize.compute_direction()
    p_new.save(path=os.path.join(tmpdir, "p_new"))

    # this query requires 'g_new' and 'p_new'
    assert(optimize.attempt_line_search_restart() == False)

    p_new.update(vector=np.array([-99.99, -99.99]))
    p_new.save(path=os.path.join(tmpdir, "p_new"))
    assert(optimize.attempt_line_search_restart() == True)


def test_line_search_recover_from_failure(tmpdir, setup_optimization_vectors):
    """
    Run a small optimization problem to try to reduce the Rosenbrock
    objective function using the Bracket'ing line search method. Simulate a job
    failure during the line search, and attempts to recover the line search
    from a checkpointed state, mimicing a real life inversion where a line
    search might fail a job (forward simulation), and we do not want to have to
    run the line search from the beginning.
    """
    optimize = Gradient(path_optimize=tmpdir, path_output=tmpdir,
                        line_search_method="bracket", step_count_max=100)

    # Make sure the initial misfit is high
    assert(optimize.load_vector("f_new") == pytest.approx(24.2, 1))

    # Calculate the initial search direction which is just the inverse gradient
    p_new = optimize.compute_direction()
    p_new.save(path=os.path.join(tmpdir, "p_new"))

    # Saves trial model 'm_try' and corresponding step length 'alpha'
    m_try, alpha = optimize.initialize_search()
    optimize.save_vector("m_try", m_try)
    optimize.save_vector("alpha", alpha)

    def line_search(_optimize, allow_break=True):
        """Each line search evaluation requires calculating a new model and
        corresponding misfit value. Simulate a break at a given step count"""
        m_try = _optimize.load_vector("m_try")
        f_try = rosenbrock_objective_function(m_try.vector)

        # Line search is most likely to break at the function evaluation
        # break mimics a job failure on cluster, before `f_try` has been saved
        if allow_break and _optimize.step_count == 3:
            return "BREAK"

        np.savetxt(os.path.join(tmpdir, "f_try.txt"), [f_try])

        m_try, alpha, status = _optimize.update_line_search()
        _optimize.save_vector("m_try", m_try)
        _optimize.save_vector("alpha", alpha)
        return status

    # Run a line search until the line search breaks
    while True:
        status = line_search(optimize, allow_break=True)
        if status == "PASS":
            break
        elif status == "BREAK":
            optimize.checkpoint()
            break

    # Try to restart the line search by re-instantiating from a checkpoint with
    # a newly instantiated optimization module
    optimize_restarted = Gradient(path_optimize=tmpdir, path_output=tmpdir,
                                  line_search_method="bracket",
                                  step_count_max=100)
    optimize_restarted.load_checkpoint()
    assert(optimize_restarted.step_count == 3)
    while True:
        status = line_search(optimize_restarted, allow_break=False)
        if status == "PASS":
            break

    # The rest of this is copied from `test_bracket_line_search` which completes
    # a successful line search. So if these values are the same then we know
    # we have successfully restarted a line search
    assert(status == "PASS")  # pass
    assert(optimize_restarted.step_count == 4)  # Took 4 steps to reduce misfit
    # Make sure we have reduced the final misfit
    assert(min(optimize_restarted._line_search.func_vals) ==
           pytest.approx(4.22, 3))

    # Change model names and reset line search
    optimize_restarted.finalize_search()

    # Check that the final model
    m_new = optimize_restarted.load_vector("m_new")
    m_old = optimize_restarted.load_vector("m_old")
    m_angle = angle(m_new.vector, m_old.vector)
    assert(m_angle == pytest.approx(0.3, 2))


def _test_inversion_optimization_problem_general(optimize, iterations=250):
    """
    Rather than run a single line search evaluation, which all the previous
    tests have done, we want to run a inversion workflow to find a best fitting
    model. To do this we essentially have to mimic the inversion workflow, but
    with barebones functions.

    This function is written to be general, other tests should populate the
    `optimize` input parameter with instantiated Optimization modules.

    .. note::
        We do not save `m_try` to disk each time it is evaluated because it is
        small. However in real workflows, `m_try` must be saved to disk rather
        than passed in memory because it is likely to be a large vector.

    .. note::
        This replaces `workflow.test_optimize` from original code

    :type optimize: module
    :param optimize: specific SeisFlows optimization module to test
    :type iterations: int
    :param iterations: number of iterations to run. defaults to 250
    """
    m_init = Model()
    m_init.model = Dict(x=[np.array([-1.2, 1])])  # Initial guess for Rosenbrock
    optimize.save_vector("m_new", m_init)

    m_true = m_init.copy()
    m_true.update(vector=np.array([1., 1.]))  # Rosenbrock global minimum

    # N allowable iterations, but reaching global min. will stop inversion
    inversion_status = "DNF"  # finished
    for iteration in range(iterations):
        # Step 1: Evaluate the objective function for given model 'm_new'
        m_new = optimize.load_vector("m_new")
        f_new = rosenbrock_objective_function(x=m_new.vector)
        optimize.save_vector("f_new", f_new)
        # Step 2: Evaluate the gradient of the objective function
        gradient = rosenbrock_gradient(x=m_new.vector)
        g_new = m_new.copy()
        g_new.update(vector=gradient)
        optimize.save_vector("g_new", g_new)
        # Step 3: Compute the search direction using the gradient
        p_new = optimize.compute_direction()
        optimize.save_vector("p_new", p_new)
        # Step 4a: Set up the line search with a trial model `m_try`
        m_try, alpha = optimize.initialize_search()
        optimize.save_vector("alpha", alpha)
        # Step 4b: Run the line search w/ various trial models 'til lower misfit
        while True:
            f_try = rosenbrock_objective_function(m_try.vector)
            optimize.save_vector("f_try", f_try)
            # Will look for the previously saved 'alpha' value
            m_try, alpha, status = optimize.update_line_search()
            if status == "PASS":
                break
            elif status == "FAIL":
                return optimize, status
            else:
                optimize.save_vector("alpha", alpha)
        optimize.save_vector("m_try", m_try)
        # Set up for the next iteration, most importantly `m_try` -> `m_new`
        optimize.finalize_search()

        # Figure out how far the updated model is from global minimum
        m_diff = np.linalg.norm(m_new.vector - m_true.vector)
        m_diff /= np.linalg.norm(m_new.vector)
        if m_diff < 1e-3:
            inversion_status = "FIN"  # finished
            break

    return optimize, inversion_status


def test_inversion_optimization_problem_with_gradient(
        tmpdir, setup_optimization_vectors):
    """Wrapper function to test the Gradient descent optimization problem"""
    gradient = Gradient(path_optimize=tmpdir, path_output=tmpdir,
                        step_count_max=40)
    optimize, status = _test_inversion_optimization_problem_general(gradient)

    # Just check a few of the stats file outputs to make sure this runs right
    assert(os.path.exists(optimize.path._stats_file))
    stats = np.genfromtxt(optimize.path._stats_file, delimiter=",", names=True)
    assert(status == "DNF")  # Did Not Finish
    assert(len(stats) == 250.)  # Fails to reach global minimum
    assert(stats["misfit"].min() == pytest.approx(0.1638, 3))
    assert(stats["if_restarted"].sum() == 0.)

    # Make plot in the tmpdir incase we need to debug the misfit reduction
    if True:
        plt.plot(stats["misfit"], "go-", markersize=2)
        plt.title("Gradient descent misfit; Rosenbrock problem")
        plt.xlabel("Iteration")
        plt.ylabel("Misfit")
        plt.axhline(1e-3, c="k")
        plt.savefig(os.path.join(tmpdir, "gradient_misfit.png"))


def test_inversion_optimization_problem_with_LBFGS(  # NOQA
        tmpdir, setup_optimization_vectors):
    """Wrapper function to test the L-BFGS descent optimization problem"""
    lbfgs = LBFGS(path_optimize=tmpdir, path_output=tmpdir,
                  line_search_method="backtrack")
    os.mkdir(lbfgs.path._LBFGS)
    optimize, status = _test_inversion_optimization_problem_general(lbfgs)

    # Just check a few of the stats file outputs to make sure this runs right
    assert(os.path.exists(optimize.path._stats_file))
    stats = np.genfromtxt(optimize.path._stats_file, delimiter=",", names=True)
    assert(status == "FIN")
    assert(len(stats) == 53.)  # reaches global min. in 95 iterations
    assert(stats["misfit"].min() == pytest.approx(1.043e-7, 3))
    assert(stats["if_restarted"].sum() == 1.)  # 1 restart

    if True:
        plt.plot(stats["misfit"], "ro-", markersize=2)
        plt.title("L-BFGS misfit; Rosenbrock problem")
        plt.xlabel("Iteration")
        plt.ylabel("Misfit")
        plt.axhline(1e-3, c="k")
        plt.savefig(os.path.join(tmpdir, "lbfgs_misfit.png"))


def test_inversion_optimization_problem_with_NLCG(  # NOQA
        tmpdir, setup_optimization_vectors):
    # NLCG will need more step counts
    nlcg = NLCG(path_optimize=tmpdir, path_output=tmpdir,
                step_count_max=40)
    optimize, status = _test_inversion_optimization_problem_general(nlcg)

    # Just check a few of the stats file outputs to make sure this runs right
    assert(os.path.exists(optimize.path._stats_file))
    stats = np.genfromtxt(optimize.path._stats_file, delimiter=",", names=True)
    assert(status == "FIN")
    assert(len(stats) == 46.)  # Fails to reach global minimum
    assert(stats["misfit"].min()) == pytest.approx(3.693E-5, 3)
    assert(stats["if_restarted"].sum() == 1.)

    if True:
        plt.plot(stats["misfit"], "bo-", markersize=2)
        plt.title("NLCG misfit; Rosenbrock problem")
        plt.xlabel("Iteration")
        plt.ylabel("Misfit")
        plt.axhline(1e-3, c="k")
        plt.savefig(os.path.join(tmpdir, "nlcg_misfit.png"))
