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

from seisflows.tools.config import Dict
from seisflows import ROOT_DIR
from seisflows.seisflows import SeisFlows
from seisflows.tools.config import load_yaml

TEST_DIR = os.path.join(ROOT_DIR, "tests")


@pytest.fixture
def par_file():
    """
    Return the test parameter file as a dictionary object
    :rtype: seisflows.config.Dict
    :return: dictionary of parameters
    """
    return os.path.join(TEST_DIR, "test_data", "parameters.yaml")


@pytest.fixture
def par_file_dict(par_file):
    """
    Return the test parameter file as a dictionary object
    :rtype: seisflows.config.Dict
    :return: dictionary of parameters
    """
    return Dict(load_yaml(par_file))


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


# def test_cmd_submit(tmpdir):
#     """
#     Test submit, also test the functionality of resume and restart which
#     are essentially wrappers for this call
#     :param tmpdir:
#     :return:
#     """
#     pass
#
#
# def test_cmd_clean(tmpdir):
#     """
#
#     :param tmpdir:
#     :return:
#     """
    pass


# def test_cmd_configure(tmpdir, setup_par_file, conf_par_file):
#     """
#     Test configuring a parameter file from a template par file
#
#     .. note::
#         I don't know exactly why, but this test needs to be run AFTER any other
#         test which runs seisflows.init(), otherwise the parameters are not
#         instantiated properly (you will hit a KeyError when trying to access
#         PAR). I think this is because of how seisflows.configure() registers
#         a relatively empty parameter file (only modules are defined), and this
#         gets saved into sys modules, affecting subsequent tests which end up
#         accessing sys.modules. I tried flushing sys.modules but it didn't work.
#         This behavior shouldn't get encountered in a real run because we
#         won't need to run init() and configure() in the same python
#         runtime environment, but I leave this warning here
#         wondering if I'll have to fix it at some point... -B
#     """
#     os.chdir(tmpdir)
#
#     # Copy in the setup par file so we can configure it
#     src = setup_par_file
#     dst = os.path.join(tmpdir, "parameters.yaml")
#     shutil.copy(src, dst)
#
#     # run seisflows init
#     with patch.object(sys, "argv", ["seisflows"]):
#         sf = SeisFlows()
#         sf.configure(relative_paths=False)
#
#     # Simple check that the configuration parameter file has the same number
#     # of lines as the one that has been created by configure
#     lines_conf = open(conf_par_file, "r").readlines()
#     lines_fill = open("parameters.yaml", "r").readlines()
#     assert (len(lines_conf) == len(lines_fill))
#
#     # My attempt to flush sys.modules which did NOT work
#     # from seisflows.tools.config import NAMES, PAR, PATH
#     # for name in NAMES:
#     #     del sys.modules[f"seisflows_{name}"]
#     # del sys.modules[PAR]
#     # del sys.modules[PATH]


def test_cmd_par(tmpdir, par_file):
    """
    Make sure the 'par' command can print and edit the parameter file
    :param tmpdir:
    :return:
    """
    # Copy given parameter file with a weird name
    src = par_file
    dst = os.path.join(tmpdir, "parameters.yaml")
    shutil.copy(src, dst)

    parameter = "workflow"
    expected_val = "forward"
    new_val = "migration"

    # testing the get option: seisflows par `parameter`
    with patch.object(sys, "argv", ["seisflows"]):
        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows", "-p", dst, "par", parameter]
        out = subprocess.run(cmd_line_arg, capture_output=True,
                             universal_newlines=True)

    # Check that we printed out the correct value
    par, val = out.stdout.strip().split(":")
    assert(par.strip() == parameter)
    assert(val.strip() == expected_val)

    # testing the set option: seisflows par `parameter` `value`
    with patch.object(sys, "argv", ["seisflows"]):
        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows","-p", dst, "par", parameter, new_val]
        out1 = subprocess.run(cmd_line_arg, capture_output=True,
                              universal_newlines=True)

        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows", "-p", dst, "par", parameter]
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

