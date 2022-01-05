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
from seisflows3.scripts.seisflows import sfparser, SeisFlows
from seisflows3.config import Dict
from seisflows3.tools.tools import loadyaml


@pytest.fixture
def par_file(tmpdir):
    """
    Copy the template parameter file into the temporary test directory
    :rtype: str
    :return: location of the parameter file
    """
    src = os.path.join("test_data", "test_parameters.yaml")
    dst = os.path.join(tmpdir, "parameters.yaml")
    shutil.copy(src, dst)
    return dst


@pytest.fixture
def par_file_dict(par_file):
    """
    Return the test parameter file as a dictionary object
    :rtype: seisflows.config.Dict
    :return: dictionary of parameters
    """
    return Dict(loadyaml(par_file))


@pytest.fixture
def sf():
    """
    Initiate an empty SeisFlows class which can be used for testing purposes
    :return:
    """

    sf = SeisFlows(return_self=True)


def test_call_seisflows(tmpdir, par_file_dict):
    """
    Make sure that calling the SeisFlows class (i.e., __init__()) works for both
    command line arguments and from within a Python environment
    """
    os.chdir(tmpdir)

    # Test calling the 'par' command from command line and inside a python
    # environemnt. Check that case-insensitivity is also honored
    # Check against the actual value coming from the parameter file
    check_val = par_file_dict.LINESEARCH
    for name in ["linesearch", "LINESEARCH", "LineSearch"]:
        # From the command line
        cmd_line_arg = ["seisflows", "par", name]
        out = subprocess.run(cmd_line_arg, capture_output=True,
                             universal_newlines=True)
        assert(out.stdout.strip() == f"{name.upper()}: {check_val}")

        # Test from inside a Python environment:
        # First we need to hackily substitute argv since our original command
        # line arguments were running the tests
        original_sys_argv = sys.argv
        sys.argv = ["python"]
        #  Then we need to redirect stdout to make sure the print statement is
        #  working as expected
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            SeisFlows(command="par", parameter=name)
        stdout = f.getvalue()
        assert(stdout.strip() == f"{name.upper()}: {check_val}")

        # Reset sys.argv just incase we want to do anything else
        sys.argv = original_sys_argv


def test_register(tmpdir, par_file_dict):
    """
    Test that the register function, which reads in PATHS and PARAMETERS
    works as expected, returning paths and parameters that we can read
    """

