"""
General test suite for all SeisFlows3 modules. Defines required parameters and
functions for each of the modules, tests importing each of the modules and that
each of the required parameters and functions exist. A sort of first-pass test
which makes sure the package is set up correctly.
"""
import os
import sys
import shutil
import pytest
from unittest.mock import patch
from seisflows3 import config
from seisflows3.seisflows import SeisFlows, return_modules


# Define dictionary dictating the bare-minimum SeisFlows3 structure.
# Each subclass will be checked to see if it meets these requirements which
# ensure that the package will work as intended. The required functions are
# determined by whether or not other submodules call for these functions, e.g.,
# an inversion workflow will call solver.eval_func()
required_structure = {
    "system": {
        "parameters": ["WALLTIME", "TASKTIME", "NTASK", "NPROC"],
        "functions": ["required", "check", "setup", "submit", "run",
                      "taskid", "checkpoint"]
    },
    "preprocess": {
        "parameters": [],
        "functions": ["required", "check", "setup", "prepare_eval_grad",
                      "sum_residuals", "finalize"]
    },
    "solver": {
        "parameters": ["MATERIALS", "DENSITY", "ATTENUATION"],
        "functions": ["required", "check", "setup", "generate_data",
                      "generate_mesh", "eval_func", "eval_grad", "load",
                      "save", "merge", "split", "source_names", "parameters"]
    },
    "postprocess": {
        "parameters": [],
        "functions": ["check", "setup", "write_gradient"]
    },
    "optimize": {
        "parameters": [],
        "functions": ["setup", "check", "compute_direction",
                      "initialize_search", "update_search",
                      "finalize_search", "retry_status",
                      "restart", "save", "load"]
    },
    "workflow": {
        "parameters": [],
        "functions": ["check", "main", "checkpoint"]
    },
}

# Just make sure that the structure here is dictated by the Config
assert(set(required_structure.keys()) == set(config.NAMES))

# Define some re-used paths
TEST_DIR = os.path.join(config.ROOT_DIR, "tests")
REPO_DIR = os.path.abspath(os.path.join(config.ROOT_DIR, ".."))


@pytest.fixture
def copy_par_file(tmpdir):
    """
    Copy the template parameter file into the temporary test directory
    :rtype: str
    :return: location of the parameter file
    """
    src = os.path.join(TEST_DIR, "test_data", "test_filled_parameters.yaml")
    dst = os.path.join(tmpdir, "parameters.yaml")
    shutil.copy(src, dst)


@pytest.fixture
def sfinit(tmpdir, copy_par_file):
    """
    Re-used function that will initate a SeisFlows3 working environment in
    sys modules
    :return:
    """
    copy_par_file
    os.chdir(tmpdir)
    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        sf._register(force=True)
    config.init_seisflows(check=False)

    return sf


def test_import(sfinit):
    """
    Test code by importing all available classes for this module.
    If any of these fails then the module itself has some code error
    (e.g., syntax errors, inheritance errors).
    """
    sfinit
    for name in config.NAMES:
        modules = return_modules()[name]
        for package, module_list in modules.items():
            for module in module_list:
                print(module)
                if "pyatoa" in module:
                    continue
                config.custom_import(name, module)()


def test_validate(sfinit, modules):
    """
    Test out path and parameter validation, essentially checking that all
    the paths and parameters are set properly

    .. note::
        This doesn't work because we have required parameters that are not
        set in the default parameter file. We can run configure beforehand
        but does that make sense?
    :return:
    """
    return
    sfinit
    for package, module_list in modules.items():
        for module in module_list:
            loaded_module = config.custom_import(MODULE, module)()
            from IPython import embed;embed()
            loaded_module.required.validate()


def test_required_parameters_exist(sfinit, modules):
    """
    Ensure that the required parameters are set in all the classes/subclasses
    That is, that the parameters defined above in REQUIRED_PARAMETERS have been
    defined by each SYSTEM class
    """
    sfinit
    for package, module_list in modules.items():
        for module in module_list:
            loaded_module = config.custom_import(MODULE, module)()
            sf_pp = loaded_module.required
            # Check that required parameters are set
            for req_par in REQUIRED_PARAMETERS:
                assert(req_par in sf_pp.parameters.keys()), \
                    f"{req_par} is a required parameter for module {MODULE}"


def test_required_functions_exist(sfinit, modules):
    """
    Make sure that the named, required functions exist within the class
    Do not execute just make sure they're defined, because they will be
    expected by other modules
    """
    sfinit
    for package, module_list in modules.items():
        for module in module_list:
            loaded_module = config.custom_import(MODULE, module)()
            for func in REQUIRED_FUNCTIONS:
                assert(func in dir(loaded_module)), \
                    f"'{func}' is a required function in module: " \
                    f"{MODULE}.{module}"


# ==============================================================================
# MODULE AND FUNCTION SPECIFIC TESTS TO FOLLOW
# ==============================================================================

