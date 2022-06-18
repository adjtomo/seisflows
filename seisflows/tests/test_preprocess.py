"""
Test suite for the SeisFlows system module, which controls interaction with
various compute systems
"""
import os
import sys
import shutil
import pytest
from unittest.mock import patch
from seisflows import config
from seisflows.seisflows import SeisFlows, return_modules


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
    Re-used function that will initate a SeisFlows working environment in
    sys modules
    :return:
    """
    # Ensure that there is not a currently active working state
    config.flush()

    copy_par_file
    os.chdir(tmpdir)
    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        sf._register_parameters(force=True)
        sf._register_modules(check=False)

    return sf


def test_default_check(sfinit):
    """
    Test seisflows.preprocess.default.check()

    :param sfinit:
    :param modules:
    :return:
    """
    sfinit
    PAR = sys.modules["seisflows_parameters"]
    preprocess = sys.modules["seisflows_preprocess"]

    # Make sure that the check statement catches incorrectly set parameters
    incorrect_parameters = {
        "NORMALIZE": ["JNORML3"],
        "MUTE": ["not_mute"],
        "FILTER": "bondpass",
        "FORMAT": "not_an_acceptable_format"
    }
    for key, val in incorrect_parameters.items():
        og_val = PAR[key]
        print(key)
        with pytest.raises(AssertionError):
            PAR[key] = val
            preprocess.check()

        PAR[key] = og_val

    # Make sure that parameters set to inappropriate values throw assertions
    correct_parameters = {
        "FILTER": "BANDPASS",
        "WORKFLOW": "INVERSION",
        "MIN_FREQ": 1,
    }
    for key, val in correct_parameters.items():
        PAR[key] = val

    incorrect_values = {
        "MAX_FREQ": -1,
    }
    for key, val in incorrect_values.items():
        og_val = PAR[key]
        with pytest.raises(AssertionError):
            PAR[key] = val
            preprocess.check()
        PAR[key] = og_val


def test_default_setup(sfinit):
    """
    Ensure that default setup correctly sets up the preprocessing machinery
    """
    sf = sfinit
    PAR = sys.modules["seisflows_parameters"]
    preprocess = sys.modules["seisflows_preprocess"]

    # Make sure that preprocess machinery is set up empty
    assert(preprocess.misfit is None)
    assert(preprocess.adjoint is None)
    assert(preprocess.reader is None)
    assert(preprocess.writer is None)

    # Set some default parameters to run setup
    misfit_name = "waveform"
    io_name = "ascii"
    PAR["MISFIT"] = misfit_name
    PAR["FORMAT"] = io_name
    preprocess.setup()

    assert(preprocess.misfit.__name__ == misfit_name)
    assert(preprocess.adjoint.__name__ == misfit_name)
    assert(preprocess.reader.__name__ == io_name)
    assert(preprocess.writer.__name__ == io_name)


def test_default_prepare_eval_grad(tmpdir, sfinit):
    """
    Ensure that prepare_eval_grad writes out adjoint traces and auxiliary files
    """
    sfinit
    PAR = sys.modules["seisflows_parameters"]
    preprocess = sys.modules["seisflows_preprocess"]

    cwd = tmpdir
    taskid = 0
    filenames = []
    preprocess.prepare_eval_grad(cwd=cwd, taskid=taskid, filenames=filenames)


