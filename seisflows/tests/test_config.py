"""
Test the SeisFlows3 configuration script, which configures the compute
system and the working environment required for SF3 to run properly
"""
import os
import sys
import shutil
import pytest
from unittest.mock import patch
from seisflows import config
from seisflows.seisflows import SeisFlows
from seisflows.tools.err import ParameterError


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
    config.init_seisflows()

    return sf


def test_seisflows_constants():
    """
    Ensure that the constants set in the Config file have not changed
    Essentially a double check to make sure these things haven't been edited
    because the rest of the package depends on these being accesible and
    the same
    """
    names_check = ["system", "preprocess", "solver",
                   "postprocess", "optimize", "workflow"]

    packages_check = ["seisflows", "seisflows-super"]

    root_dir_check = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), ".."
    )

    assert(config.NAMES == names_check)
    assert(config.PACKAGES == packages_check)
    assert(os.path.samefile(config.ROOT_DIR, root_dir_check))


def test_init_seisflows(sfinit):
    """
    Make sure that initiation of the modular approach of seisflows works
    as expected. That is, that system-wide accessible modules are
    instantiated with the accepted naming schema

    .. note::
        This assumes that the parameter file is set up correctly, which it
        should be if it's coming from the test data directory
    :return:
    """
    sf = sfinit
    # Ensure that all the modules in NAMES have been instantiated in sys.modules
    for name in config.NAMES:
        assert(f"seisflows_{name}" in sys.modules)


def test_save_and_load(sfinit):
    """
    Test saving the current session to disk
    :return:
    """
    # Instantiate sys modules and save to disk
    sf = sfinit
    config.save()
    # Now remove seisflows sys modules so we can try load them back
    for name in config.NAMES:
        sys.modules.pop(f"seisflows_{name}")
    config.load(path="./output")
    for name in config.NAMES:
        assert(f"seisflows_{name}" in sys.modules)


def test_seisflows_paths_parameters(sfinit):
    """
    Test the class that makes inputting and checking paths and parameters easier
    Recreates the required() function at the top of each class.
    """
    sf = sfinit
    sfpp = config.SeisFlowsPathsParameters()

    # All of these parameters are defined in the test parameter file
    sfpp.par("SOLVER", required=True, par_type=str,
           docstr="This is a required parameter")
    sfpp.par("MIN_PERIOD", required=False, default=10., par_type=float,
           docstr="This is an optional parameter")
    sfpp.path("SPECFEM_BIN", required=True, docstr="This is a required path")
    sfpp.path("LOCAL", required=False, docstr="This is an optional path")
    sfpp.validate()

    # These parameters are not defined and are expected to throw parameter error
    sfpp.path("UNDEFINED", required=True,
              docstr="This path is not in the test parameter file")
    with pytest.raises(ParameterError):
        sfpp.validate()


def test_custom_import(sfinit):
    """
    Test that importing based on internal modules works for various inputs
    :return:
    """
    sfinit
    with pytest.raises(SystemExit):
        config.custom_import()
    with pytest.raises(SystemExit):
        config.custom_import(name="NOT A VALID NAME")

    module = config.custom_import(name="optimize", module="LBFGS")
    assert(module.__name__ == "LBFGS")
    assert(module.__module__ == "seisflows.optimize.LBFGS")

    # Check one more to be safe
    module = config.custom_import(name="optimize", module="base")
    assert(module.__name__ == "Base")
    assert(module.__module__ == "seisflows.optimize.base")


