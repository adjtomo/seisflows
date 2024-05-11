"""
Test any of the utility functions defined in the Tools directory
"""
import os
import pytest
import numpy as np
from glob import glob
from seisflows import ROOT_DIR
from seisflows.tools.config import Dict
from seisflows.tools.model import Model
from seisflows.tools.config import custom_import


TEST_DIR = os.path.join(ROOT_DIR, "tests")
TEST_MODEL = os.path.join(TEST_DIR, "test_data", "test_tools", 
                          "test_file_formats")


def test_model_attributes():
    """
    Check that model values are read in correctly and accessible in a way we 
    expect
    """
    m = Model(path=TEST_MODEL, fmt=".bin")

    assert(m.ngll[0] == 40000)
    assert(m.nproc == 1)
    assert("vp" in m.model.keys())
    assert("vs" in m.model.keys())
    assert(m.model.vp[0][0] == 5800.)
    assert(m.model.vs[0][0] == 3500.)


def test_model_functions():
    """
    Test the core merge and split functions of the Model class to ensure
    they give the same results
    """
    m = Model(path=TEST_MODEL, fmt=".bin")

    assert(len(m.merge() == len(m.model.vs[0]) + len(m.model.vp[0])))
    assert(len(m.split()) == len(m.parameters))


def test_model_io(tmpdir):
    """
    Check that saving and loading npz file works
    """
    m = Model(path=TEST_MODEL, fmt=".bin")

    m.save(path=os.path.join(tmpdir, "test.npz"))
    m_new = Model(path=os.path.join(tmpdir, "test.npz"))
    assert(m_new.ngll[0] == m.ngll[0])
    assert(m_new.fmt == m.fmt)

    # Check that writing fortran binary works
    m.write(path=tmpdir)
    assert(len(glob(os.path.join(tmpdir, f"*{m.fmt}"))) == 2)


def test_model_from_input_vector():
    """Check that we can instantiate a model from an input vector"""
    m = Model(path=None)
    m.model = Dict(x=[np.array([-1.2, 1.])])
    assert(m.nproc == 1)
    assert(m.ngll == [2])
    assert(m.parameters == ["x"])


def test_custom_import():
    """
    Test that importing based on internal modules works for various inputs
    :return:
    """
    with pytest.raises(SystemExit):
        custom_import()
    with pytest.raises(SystemExit):
        custom_import(name="NOT A VALID NAME")

    module = custom_import(name="optimize", module="LBFGS")
    assert(module.__name__ == "LBFGS")
    assert(module.__module__ == "seisflows.optimize.LBFGS")

    # Check one more to be safe
    module = custom_import(name="preprocess", module="default")
    assert(module.__name__ == "Default")
    assert(module.__module__ == "seisflows.preprocess.default")
