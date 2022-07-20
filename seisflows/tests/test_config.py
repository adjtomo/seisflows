"""
Test the SeisFlows configuration script, which configures the compute
system and the working environment required for SF to run properly
"""
import os
import shutil
import pytest

from seisflows import config


TEST_DIR = os.path.join(config.ROOT_DIR, "tests")


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


def test_seisflows_constants():
    """
    Ensure that the constants set in the Config file have not changed
    Essentially a double check to make sure these things haven't been edited
    because the rest of the package depends on these being accesible and
    the same
    """
    names_check = ["system", "preprocess", "solver", "optimize", "workflow"]

    root_dir_check = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), ".."
    )

    assert(config.NAMES == names_check)
    assert(os.path.samefile(config.ROOT_DIR, root_dir_check))


def test_custom_import():
    """
    Test that importing based on internal modules works for various inputs
    :return:
    """
    with pytest.raises(SystemExit):
        config.custom_import()
    with pytest.raises(SystemExit):
        config.custom_import(name="NOT A VALID NAME")

    module = config.custom_import(name="optimize", module="LBFGS")
    assert(module.__name__ == "LBFGS")
    assert(module.__module__ == "seisflows.optimize.LBFGS")

    # Check one more to be safe
    module = config.custom_import(name="preprocess", module="default")
    assert(module.__name__ == "Default")
    assert(module.__module__ == "seisflows.preprocess.default")


