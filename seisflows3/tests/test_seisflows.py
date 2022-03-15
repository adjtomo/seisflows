"""
Test suite for the SeisFlows command line interface tool and the underlying
parser class, ensures that the command line tool works as expected

!!! TO DO: Finish tests from 'edit' onwards
"""
import os
import io
import sys
import shutil
import pytest
import contextlib
import subprocess
from unittest.mock import patch
from seisflows3.seisflows import sfparser, SeisFlows
from seisflows3.config import Dict, ROOT_DIR
from seisflows3.tools.wrappers import loadyaml

TEST_DIR = os.path.join(ROOT_DIR, "tests")


@pytest.fixture
def filled_par_file():
    """A parameter file that is completely filled out and can be read in"""
    return os.path.join(TEST_DIR, "test_data", "test_filled_parameters.yaml")


@pytest.fixture
def conf_par_file():
    """A par file that has been configured but requires user-defined values"""
    return os.path.join(TEST_DIR, "test_data", "test_conf_parameters.yaml")


@pytest.fixture
def setup_par_file():
    """A barebones par file that only contains module names"""
    return os.path.join(TEST_DIR, "test_data", "test_setup_parameters.yaml")


@pytest.fixture
def copy_par_file(tmpdir, filled_par_file):
    """
    Copy the template parameter file into the temporary test directory
    :rtype: str
    :return: location of the parameter file
    """
    src = filled_par_file
    dst = os.path.join(tmpdir, "parameters.yaml")
    shutil.copy(src, dst)


@pytest.fixture
def par_file_dict(filled_par_file):
    """
    Return the test parameter file as a dictionary object
    :rtype: seisflows.config.Dict
    :return: dictionary of parameters
    """
    return Dict(loadyaml(filled_par_file))


def test_call_seisflows(tmpdir, par_file_dict, copy_par_file):
    """
    Test calling the 'par' command from command line and inside a python
    environemnt. Check that case-insensitivity is also honored
    Check against the actual value coming from the parameter file
    Also tests the seisflows 'par' command
    """
    copy_par_file
    os.chdir(tmpdir)
    check_val = par_file_dict.LINESEARCH

    # Mock argv to match the actual scenario in which SeisFlows will be called
    with patch.object(sys, "argv", ["seisflows"]):
        for name in ["linesearch", "LINESEARCH", "LineSearch"]:
            # From the command line
            cmd_line_arg = ["seisflows", "par", name]
            out = subprocess.run(cmd_line_arg, capture_output=True,
                                 universal_newlines=True)
            assert(out.stdout.strip() == f"{name.upper()}: {check_val}")

            # Test from inside a Python environment; we need to redirect stdout to
            # make sure the print statement is working as expected
            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                sf = SeisFlows()
                sf(command="par", parameter=name)
            stdout = f.getvalue()
            assert(stdout.strip() == f"{name.upper()}: {check_val}")


def test_edited_parameter_file_name(tmpdir, par_file_dict, filled_par_file):
    """
    Similar test as call_seisflows but just make sure that arbitrary naming
    of the parameter file still works
    """
    par_fid = "CRAZY_PARAMETER_NAME-0x123kd.yaml"
    par_name = "linesearch"

    src = filled_par_file
    dst = os.path.join(tmpdir, par_fid)
    shutil.copy(src, dst)

    os.chdir(tmpdir)
    check_val = par_file_dict.LINESEARCH
    cmd_line_arg = ["seisflows", "-p", par_fid, "par", par_name]
    # Mock argv to match the actual scenario in which SeisFlows will be called
    with patch.object(sys, "argv", ["seisflows"]):
        out = subprocess.run(cmd_line_arg, capture_output=True,
                             universal_newlines=True)
    assert(out.stdout.strip() == f"{par_name.upper()}: {check_val}")


def test_register(tmpdir, par_file_dict, copy_par_file):
    """
    Test that the register function, which reads in PATHS and PARAMETERS
    works as expected, returning paths and parameters that we can read
    """
    copy_par_file
    os.chdir(tmpdir)

    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        assert(sf._paths is None)
        assert(sf._parameters is None)
        sf._register(force=True)

    # Check that paths and parameters have been set in sys.modules
    paths = sys.modules["seisflows_paths"]
    parameters = sys.modules["seisflows_parameters"]

    # Check one or two parameters have been set correctly
    assert(par_file_dict.LBFGSMAX == parameters.LBFGSMAX)
    path_check_full = os.path.abspath(par_file_dict.PATHS["SCRATCH"])
    assert(path_check_full == paths.SCRATCH)


def test_cmd_setup(tmpdir):
    """
    Test setting up the SeisFlows3 working directory
    """
    os.chdir(tmpdir)
    par_file = os.path.join(tmpdir, "parameters.yaml")

    # Check with anmd without symlinking as well as with overwriting
    with patch.object(sys, "argv", ["seisflows"]):
        # Without symlinking
        sf = SeisFlows()
        sf.setup(symlink=False, overwrite=False)
        assert(os.path.exists(par_file))
        os.remove(par_file)

        # With symlinking
        sf.setup(symlink=True, overwrite=False)
        assert(os.path.exists(par_file))
        assert(os.path.exists(
            os.path.join(tmpdir, "source_code", "seisflows3"))
        )
        # Edit the current par file in a noticeable way so we can check
        # if overwriting works in the next step
        test_phrase = "well this is rather unexpected...\n"
        with open(par_file, "a") as f:
            f.write(test_phrase)
        with open(par_file, "r") as f:
            assert(test_phrase in f.read())

        # With overwriting
        sf.setup(symlink=False, overwrite=True)
        assert(os.path.exists(par_file))
        with open(par_file, "r") as f:
            text = f.read()
            assert(test_phrase not in f.read())


def test_cmd_configure(tmpdir, setup_par_file, conf_par_file):
    """
    Test configuring a parameter file
    """
    os.chdir(tmpdir)

    # Copy in the setup par file so we can configure it
    src = setup_par_file
    dst = os.path.join(tmpdir, "parameters.yaml")
    shutil.copy(src, dst)

    # Configure the empty parameter file
    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        # Need to run this via command line because there are sub arguments
        cmd_line_arg = ["seisflows", "configure"]
        subprocess.run(cmd_line_arg, capture_output=True,
                       universal_newlines=True)

    # We will be checking against an already configured par file
    for check, par in zip(open(conf_par_file).readlines(),
                          open(dst).readlines()):
        # Titles will be difference since it's based on relative path
        if "TITLE" in check:
            continue
        # Paths will also be different since they're relative
        elif "PATH" in check:
            break
        # Otherwise all the lines should be the same
        assert(check == par)


def blank(tmpdir):
    """
    Test setting up the SeisFlows3 working directory
    """
    pass
    os.chdir(tmpdir)

    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        sf._register(precheck=False)

