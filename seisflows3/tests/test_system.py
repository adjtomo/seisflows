"""
Test suite for the SeisFlows3 system module, which controls interaction with
various compute systems
"""
import os
import sys
import shutil
import pytest
from unittest.mock import patch
from seisflows3 import config
from seisflows3.seisflows import SeisFlows, return_modules


# The module that we're testing, allows for copy-pasting these test suites
MODULE = "system"

# Ensures that these parameters are always defined, even when using subclasses
REQUIRED_PARAMETERS = ["WALLTIME", "TASKTIME", "NTASK", "NPROC"]
REQUIRED_FUNCTIONS = ["required", "check", "setup", "submit", "run",
                      "run_single", "taskid", "checkpoint", "wrong"
                      ]

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
def modules():
    """
    Return a list of subclasses that falls under the System module
    """
    return return_modules()[MODULE]


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


def test_import(sfinit, modules):
    """
    Test code by importing all available classes for this module.
    If any of these fails then the module itself has some code error
    (e.g., syntax errors, inheritance errors).
    """
    sfinit
    for package, module_list in modules.items():
        for module in module_list:
            config.custom_import(MODULE, module)()


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
    pass
    sfinit
    for package, module_list in modules.items():
        for module in module_list:
            loaded_module = config.custom_import(MODULE, module)()
            pytest.set_trace()
            from IPython import embed;embed()
            loaded_module.required.validate()



def test_required_parameters(sfinit, modules):
    """
    Ensure that the required parameters are set in all the classes/subclasses
    """
    sf = sfinit
    for package, module_list in modules.items():
        for module in module_list:
            loaded_module = config.custom_import(MODULE, module)()
            sf_pp = loaded_module.required
            # Check that required parameters are set
            for req_par in REQUIRED_PARAMETERS:
                assert(req_par in sf_pp.parameters.keys()), \
                    f"{req_par} is a required parameter for module {MODULE}"

            pytest.set_trace()


def test_required_functions(sfinit, modules):
    """
    Make sure that the named, required functions exist within the class
    """
    sf = sfinit
    for package, module_list in modules.items():
        for module in module_list:
            loaded_module = config.custom_import(MODULE, module)()
            for func in REQUIRED_FUNCTIONS:
                assert(func in dir(loaded_module)), \
                    f"'{func}' is a required function in module: {MODULE}.{module}"


def test_setup(sfinit, modules):
    """
    Make sure that setup creates the necessary directory structure

    :param sfinit:
    :param modules:
    :return:
    """
    sf = sfinit
    PATHS = sys.modules["seisflows_paths"]
    SETUP_CREATES = [PATHS.SCRATCH, PATHS.SYSTEM, PATHS.OUTPUT]
    for package, module_list in modules.items():
        for module in module_list:
            loaded_module = config.custom_import(MODULE, module)()

            # Make sure these don't already exist
            for path_ in SETUP_CREATES:
                assert(not os.path.exists(path_))

            loaded_module.setup()

            # Check that the minimum required directories were created
            for path_ in SETUP_CREATES:
                assert(os.path.exists(path_))

            # Remove created paths so we can check the next module
            for path_ in SETUP_CREATES:
                if os.path.isdir(path_):
                    shutil.rmtree(path_)
                else:
                    os.remove(path_)


def test_submit():
    """

    :return:
    """


def test_run():
    """

    :return:
    """


def test_run_single():
    """

    :return:
    """


def test_taskid():
    """
    Simply assert that this function returns an integer
    :return:
    """

