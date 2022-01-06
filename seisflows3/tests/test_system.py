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
from seisflows3.scripts.seisflows import SeisFlows, return_modules


# The module that we're testing, allows for copy-pasting these test suites
MODULE = "system"

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
        sf._register(precheck=False)
    config.init_seisflows()

    return sf


def test_import(sfinit, modules):
    """
    Test code execution by importing all of the available modules in the package
    If any of these fails then the module itself has some error (e.g.,
    syntax errors) or the 'required' statement is failing
    """
    sf = sfinit
    for package, module_list in modules.items():
        for module in module_list:
            loaded_module = config.custom_import(MODULE, module)()
            # Not sure what the best way to check these things are but
            # for now we run a validate which just makes sure all the
            # paths and parameters are set into sys.modules
            loaded_module.required.validate()


def test_required(sfinit, modules):
    """
    Ensure that the path/parameter checking mechanism is working as advertised
    """
    REQUIRED_PARAMETERS = ["WALLTIME", "TASKTIME", "NTASK", "NPROC"]
    sf = sfinit

    pytest.set_trace()
