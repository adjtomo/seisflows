"""
Test suite for the SeisFlows3 SYSTEM module, which controls interaction with
various compute systems
"""
import os
import sys
import shutil
import pytest
from glob import glob
from unittest.mock import patch
from seisflows import config
from seisflows.seisflows import SeisFlows


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


def test_setup(sfinit):
    """
    Make sure that system.base.setup() creates the desired directory structure
    and log files.

    :param sfinit:
    :param modules:
    :return:
    """
    sf = sfinit
    system = sys.modules["seisflows_system"]
    PATH = sys.modules["seisflows_paths"]

    # Make empty text files to check that they will be copied into the logs dir
    for fid in [system.output_log, system.error_log]:
        open(fid, "w")

    for path in [PATH.SCRATCH, PATH.SYSTEM, PATH.OUTPUT]:
        assert(not os.path.exists(path))

    system.setup()
    for path in [PATH.SCRATCH, PATH.SYSTEM, PATH.OUTPUT]:
        assert(os.path.exists(path))

    # Both log files and the parameter file should have been copied
    assert(len(glob(os.path.join("logs", "*"))) == 3)

def test_checkpoint(sfinit):
    """
    Check that output pickle files are created during a checkpoint and that
    kwargs are saved
    """
    sf = sfinit
    system = sys.modules["seisflows_system"]
    PATH = sys.modules["seisflows_paths"]
    classname = "solver"
    method = "eval_func"
    system.checkpoint(path=PATH.OUTPUT, classname=classname, method=method,
                      kwargs={"test": 5}
                      )
    assert(os.path.exists(os.path.join(PATH.OUTPUT, "kwargs",
                                       f"{classname}_{method}.p")))
    assert(len(glob(os.path.join(PATH.OUTPUT, "*.p"))) == len(config.NAMES))


# SYSTEM.WORKSTATION
def test_run_system_workstation(sfinit):
    """
    Ensure that workstation.run() simply runs a function as we would expect
    """
    sf = sfinit
    system = sys.modules["seisflows_system"]
    assert(type(system).__name__ == "Workstation")
    # We don't care what the function does, just that we can call it
    system.run(classname="system", method="taskid")

