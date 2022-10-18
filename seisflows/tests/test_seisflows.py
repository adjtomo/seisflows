"""
Test suite for the SeisFlows command line interface tool and the underlying
parser class, ensures that the command line tool works as expected
"""
import os
import io
import sys
import shutil
import pytest
import contextlib
import subprocess
from unittest.mock import patch

from seisflows import ROOT_DIR
from seisflows.tools import msg
from seisflows.seisflows import SeisFlows
from seisflows.tools.config import load_yaml, Dict

TEST_DIR = os.path.join(ROOT_DIR, "tests")


@pytest.fixture
def par_file(tmpdir):
    """
    Create a template parameter file
    """
    fid = os.path.join(tmpdir, "parameters.yaml")
    with open(fid, "w") as f:
        f.write(msg.base_parameter_file)
    return fid


@pytest.fixture
def par_file_dict():
    """
    Return default parameter file parameters as a dictionary object

    :rtype: seisflows.config.Dict
    :return: dictionary of parameters
    """
    par_file = Dict(workflow="forward", system="workstation",
                    solver="specfem2d", preprocess="default",
                    optimize="gradient")
    return par_file


def test_call_seisflows(par_file, par_file_dict):
    """
    Test calling the 'par' command from command line and inside a python
    environemnt. Check that case-insensitivity is also honored
    Check against the actual value coming from the parameter file
    Also tests the seisflows 'par' command
    """
    check_val = par_file_dict.workflow

    # Mock argv to match the actual scenario in which SeisFlows will be called
    with patch.object(sys, "argv", ["seisflows"]):
        for name in ["Workflow", "WORKFLOW", "WorkFlow"]:
            # $ seisflows par workflow -p path/to/parameters.yaml
            cmd_line_arg = ["seisflows", "-p", par_file, "par", name]
            out = subprocess.run(cmd_line_arg, capture_output=True,
                                 universal_newlines=True)
            assert(out.stdout.strip() == f"{name.lower()}: {check_val}")

            # Test from inside a Python environment; we need to redirect stdout
            # to make sure the print statement is working as expected
            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                sf = SeisFlows(parameter_file=par_file)
                sf(command="par", parameter=name)
            stdout = f.getvalue()
            assert(stdout.strip() == f"{name.lower()}: {check_val}")


def test_edited_parameter_file_name(tmpdir, par_file, par_file_dict):
    """
    Similar test as call_seisflows but just make sure that arbitrary naming
    of the parameter file still works
    """
    # Copy given parameter file with a weird name
    par_fid = "CRAZY_PARAMETER_NAME-0x123kd.yaml"
    src = par_file
    dst = os.path.join(tmpdir, par_fid)
    shutil.copy(src, dst)
    check_key = "workflow"
    check_val = par_file_dict[check_key]

    cmd_line_arg = ["seisflows", "-p", dst, "par", check_key]
    # Mock argv to match the actual scenario in which SeisFlows will be called
    with patch.object(sys, "argv", ["seisflows"]):
        out = subprocess.run(cmd_line_arg, capture_output=True,
                             universal_newlines=True)
    assert(out.stdout.strip() == f"{check_key.lower()}: {check_val}")


def test_cmd_setup(tmpdir):
    """
    Test setting up the SeisFlows working directory
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

        # Edit the current par file in a noticeable way so we can check
        # if overwriting works in the next step
        test_phrase = "well this is rather unexpected...\n"
        with open(par_file, "a") as f:
            f.write(test_phrase)
        with open(par_file, "r") as f:
            assert(test_phrase in f.read())

        # With overwriting
        sf.setup(symlink=False, overwrite=True, force=True)
        assert(os.path.exists(par_file))
        with open(par_file, "r") as f:
            text = f.read()
            assert(test_phrase not in text)


def test_cmd_configure(tmpdir, par_file):
    """
    Test configuring a parameter file from a template par file
    """
    os.chdir(tmpdir)

    # run seisflows configure
    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows(workdir=tmpdir, parameter_file="parameters.yaml")
        sf.configure()

    # Check some random values that were not in the template file
    parameters = load_yaml(par_file)
    assert("path_model_init" in parameters.keys())
    assert("smooth_h" in parameters.keys())
    assert("ntask" in parameters.keys())


def test_cmd_par(tmpdir, par_file):
    """
    Make sure the 'par' command can print and edit the parameter file
    :param tmpdir:
    :return:
    """
    parameter = "workflow"
    expected_val = "forward"
    new_val = "migration"

    # testing the get option: seisflows par `parameter`
    with patch.object(sys, "argv", ["seisflows"]):
        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows", "-p", par_file, "par", parameter]
        out = subprocess.run(cmd_line_arg, capture_output=True,
                             universal_newlines=True)

    # Check that we printed out the correct value
    par, val = out.stdout.strip().split(":")
    assert(par.strip() == parameter)
    assert(val.strip() == expected_val)

    # testing the set option: seisflows par `parameter` `value`
    with patch.object(sys, "argv", ["seisflows"]):
        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows", "-p", par_file, "par", parameter, new_val]
        out1 = subprocess.run(cmd_line_arg, capture_output=True,
                              universal_newlines=True)

        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows", "-p", par_file, "par", parameter]
        out2 = subprocess.run(cmd_line_arg, capture_output=True,
                              universal_newlines=True)

    # Check that the changed print statement works
    par, vals = out1.stdout.strip().split(":")
    val_old, val_new = vals.strip().split(" -> ")
    assert(par.strip() == parameter)
    assert(val_old.strip() == expected_val)
    assert(val_new.strip() == new_val)

    # Check that we printed out the correctly changed value
    par, val = out2.stdout.strip().split(":")
    assert(par.strip() == parameter)
    assert(val.strip() == new_val)

