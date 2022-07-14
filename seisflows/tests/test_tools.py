"""
Test any of the utility functions defined in the Tools directory
"""
import os
import pytest
from glob import glob
from seisflows.tools.specfem import Model
from seisflows.config import ROOT_DIR, NAMES, CFGPATHS


TEST_DIR = os.path.join(ROOT_DIR, "tests")


def test_specfem_model(tmpdir):
    """
    Make sure we can dynamically load SPECFEM models in various formats

    TODO eventually we want to split this into multiple tests that test each
    TODO of the IO formats (binary, adios etc.). Currently only testing binary.
    """
    model_data = os.path.join(TEST_DIR, "test_data", "test_tools",
                              "test_file_formats")

    # Make sure that multiple acceptable file extensions throw error
    with pytest.raises(AssertionError):
        Model(path=model_data)

    # Check that model values are read in correctly
    m = Model(path=model_data, fmt=".bin")
    assert(m.ngll[0] == 40000)
    assert(m.nproc == 1)
    assert("vp" in m.model.keys())
    assert("vs" in m.model.keys())
    assert(m.model.vp[0][0] == 5800.)
    assert(m.model.vs[0][0] == 3500.)

    assert(len(m.merge() == len(m.model.vs[0]) + len(m.model.vp[0])))
    assert(len(m.split()) == len(m.parameters))

    # Check that saving and loading npz file works
    m.save(path=os.path.join(tmpdir, "test.npz"))
    m_new = Model(path=os.path.join(tmpdir, "test.npz"), load=True)
    assert(m_new.ngll[0] == m.ngll[0])
    assert(m_new.fmt == m.fmt)

    # Check that writing fortran binary works
    m.write(path=tmpdir)
    assert(len(glob(os.path.join(tmpdir, f"*{m.fmt}"))) == 2)
