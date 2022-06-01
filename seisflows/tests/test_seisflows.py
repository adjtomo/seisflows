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

from seisflows import logger
from seisflows.seisflows import sfparser, SeisFlows
from seisflows.config import (save, Dict, ROOT_DIR, NAMES, CFGPATHS,
                              config_logger)
from seisflows.tools.wrappers import loadyaml

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

            # Test from inside a Python environment; we need to redirect stdout
            # to make sure the print statement is working as expected
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
        assert(os.path.exists(
            os.path.join(tmpdir, "source_code", "seisflows"))
        )
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


def test_cmd_init(tmpdir, copy_par_file):
    """
    Test 'seisflows init' command which instantiates a working directory and
    saves the active working state as pickle files
    :return:
    """
    os.chdir(tmpdir)
    copy_par_file

    # Create necessary paths to get past some assertion errors
    parameters = loadyaml("parameters.yaml")
    paths = parameters.pop("PATHS")

    for key in ["MODEL_INIT", "MODEL_TRUE"]:
        os.mkdir(paths[key])

    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        sf.init()

    for name in NAMES:
        assert(os.path.exists(os.path.join(paths["OUTPUT"],
                                           f"seisflows_{name}.p"))
               )


def test_cmd_submit(tmpdir):
    """
    Test submit, also test the functionality of resume and restart which
    are essentially wrappers for this call
    :param tmpdir:
    :return:
    """
    pass


def test_cmd_clean(tmpdir):
    """

    :param tmpdir:
    :return:
    """
    os.chdir(tmpdir)

    # Create a bunch of files that match what should be deleted. Make them as
    # directories even though some should be files because we just want to see
    # if they get deleted or not
    for path in CFGPATHS.values():
        os.mkdir(path)

    # Symlink the last file to make sure it still exists even if it matches
    shutil.rmtree(path)
    os.symlink(src=CFGPATHS.PAR_FILE, dst=path)

    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        sf.clean(force=True)

    for fid in [path, CFGPATHS.PAR_FILE]:
        assert(os.path.exists(fid))


def test_config_logging(tmpdir, copy_par_file):
    """
    Test logging configuration to make sure we can print to file

    TODO move this to test_config.py?
    :param tmpdir:
    :return:
    """
    # Run init first to create a working state
    os.chdir(tmpdir)
    copy_par_file

    msg = "This is an example log that will be checked for test purposes"
    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        sf._register(force=True)
        config_logger(filename=CFGPATHS.LOGFILE)
        logger.debug(msg)

    # Check that we created the log file and wrote the message in
    assert(os.path.exists(CFGPATHS.LOGFILE))
    with open(CFGPATHS.LOGFILE, "r") as f:
        lines = f.read()
    assert(msg in lines)
    assert("DEBUG" in lines)  # levelname
    assert("test_config_logging()" in lines)  # funcName


def test_load_modules(tmpdir, copy_par_file):
    """
    Test if module loading from sys.modules works

    :param tmpdir:
    :return:
    """
    # Run init first to create a working state
    os.chdir(tmpdir)
    copy_par_file

    # Create necessary paths to get past some assertion errors
    parameters = loadyaml("parameters.yaml")
    paths = parameters.pop("PATHS")

    for key in ["MODEL_INIT", "MODEL_TRUE"]:
        os.mkdir(paths[key])

    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        sf.init()

    # Check a random parameter and then set it to something different
    preprocess = sys.modules["seisflows_preprocess"]
    assert(preprocess.misfit is None)
    preprocess.misfit = 1

    # See if we can load modules and restore previous working state which
    # overwrites the previous operation
    sf._load_modules()
    assert(sys.modules["seisflows_preprocess"].misfit != 1)


def test_cmd_configure(tmpdir, setup_par_file, conf_par_file):
    """
    Test configuring a parameter file from a template par file

    .. note::
        I don't know exactly why, but this test needs to be run AFTER any other
        test which runs seisflows.init(), otherwise the parameters are not
        instantiated properly (you will hit a KeyError when trying to access
        PAR). I think this is because of how seisflows.configure() registers
        a relatively empty parameter file (only modules are defined), and this
        gets saved into sys modules, affecting subsequent tests which end up
        accessing sys.modules. I tried flushing sys.modules but it didn't work.
        This behavior shouldn't get encountered in a real run because we
        won't need to run init() and configure() in the same python
        runtime environment, but I leave this warning here
        wondering if I'll have to fix it at some point... -B
    """
    os.chdir(tmpdir)

    # Copy in the setup par file so we can configure it
    src = setup_par_file
    dst = os.path.join(tmpdir, "parameters.yaml")
    shutil.copy(src, dst)

    # run seisflows init
    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        sf.configure(relative_paths=False)

    # Simple check that the configuration parameter file has the same number
    # of lines as the one that has been created by configure
    lines_conf = open(conf_par_file, "r").readlines()
    lines_fill = open("parameters.yaml", "r").readlines()
    assert (len(lines_conf) == len(lines_fill))

    # My attempt to flush sys.modules which did NOT work
    # from seisflows.config import NAMES, PAR, PATH
    # for name in NAMES:
    #     del sys.modules[f"seisflows_{name}"]
    # del sys.modules[PAR]
    # del sys.modules[PATH]


def test_cmd_par(tmpdir, copy_par_file):
    """
    Make sure the 'par' command can print and edit the parameter file
    :param tmpdir:
    :return:
    """
    # Run init first to create a working state
    os.chdir(tmpdir)
    copy_par_file

    parameter = "begin"
    expected_val = "1"
    new_val = "2"

    # testing the get option: seisflows par `parameter`
    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows", "par", parameter]
        out = subprocess.run(cmd_line_arg, capture_output=True,
                             universal_newlines=True)

    # Check that we printed out the correct value
    par, val = out.stdout.strip().split(":")
    assert(par.upper() == parameter.upper())
    assert(int(val) == int(expected_val))

    # testing the set option: seisflows par `parameter` `value`
    with patch.object(sys, "argv", ["seisflows"]):
        sf = SeisFlows()
        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows", "par", parameter, new_val]
        out1 = subprocess.run(cmd_line_arg, capture_output=True,
                              universal_newlines=True)

        # Run this with subprocess so we can capture the print statement
        cmd_line_arg = ["seisflows", "par", parameter]
        out2 = subprocess.run(cmd_line_arg, capture_output=True,
                              universal_newlines=True)

    # Check that the changed print statement works
    par, vals = out1.stdout.strip().split(":")
    val_old, val_new = vals.strip().split(" -> ")
    assert(par.upper() == parameter.upper())
    assert(int(val_old) == int(expected_val))
    assert(int(val_new) == int(new_val))

    # Check that we printed out the correctly changed value
    par, val = out2.stdout.strip().split(":")
    assert(par.upper() == parameter.upper())
    assert(int(val) == int(new_val))

# def test_cmd_sempar(tmpdir):
#     """
#
#     :param tmpdir:
#     :return:
#     """
#     pass
#
#
# def test_cmd_check(tmpdir):
#     """
#     Very simple
#     :param tmpdir:
#     :return:
#     """
#     pass
#
#
# def test_cmd_print(tmpdir):
#     """
#
#     :param tmpdir:
#     :return:
#     """
#     pass
#
#
# def test_cmd_convert(tmpdir):
#     """
#
#     :param tmpdir:
#     :return:
#     """
#     pass
#
#
# def test_cmd_validate(tmpdir):
#     """
#
#     :param tmpdir:
#     :return:
#     """
#     pass
