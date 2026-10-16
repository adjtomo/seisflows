"""
Test any of the utility functions defined in the Tools directory
"""
import os
import pytest
import numpy as np
from glob import glob
from obspy import Trace
from seisflows import ROOT_DIR
from seisflows.tools.config import Dict
from seisflows.tools.model import Model
from seisflows.tools.config import custom_import, atomic_write_lines
from seisflows.tools.noise import rotate_ne_trace_to_rt, rotate_rt_adjsrc_to_ne


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
    Check that writing and re-reading a model in Fortran binary format works
    (round-trips through disk the same way SeisFlows <-> SPECFEM does)
    """
    m = Model(path=TEST_MODEL, fmt=".bin")

    # `Model` only knows how to write real SPECFEM-recognized parameter names
    # in Fortran binary format, so round-trip through that rather than a
    # (no-longer-supported) generic npz save
    m.write(path=tmpdir)
    assert(len(glob(os.path.join(tmpdir, f"*{m.fmt}"))) == 2)

    m_new = Model(path=tmpdir, fmt=".bin")
    assert(m_new.ngll[0] == m.ngll[0])
    assert(m_new.fmt == m.fmt)
    assert(np.asarray(m_new.vector, dtype=float) ==
           pytest.approx(np.asarray(m.vector, dtype=float), rel=1e-6))


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


def _trace(data):
    """Convenience function to build a minimal ObsPy Trace for rotation tests"""
    return Trace(data=np.array(data, dtype=float))


def test_rotate_ne_trace_to_rt_theta_zero():
    """
    At theta==theta_p==0, Wang et al. (2019) Eqs. 9/10 reduce to TT==EE and
    RR==NN (all other terms drop out via sin(0)==0/cos(0)==1)
    """
    tr_ee, tr_ne = _trace([1., 2., 3.]), _trace([10., 20., 30.])
    tr_en, tr_nn = _trace([100., 200., 300.]), _trace([4., 5., 6.])

    tr_rr, tr_tt = rotate_ne_trace_to_rt(tr_ee=tr_ee, tr_ne=tr_ne, tr_en=tr_en,
                                         tr_nn=tr_nn, theta=0., theta_p=0.)

    assert(tr_tt.data == pytest.approx(tr_ee.data))
    assert(tr_rr.data == pytest.approx(tr_nn.data))
    assert(tr_tt.stats.component == "T")
    assert(tr_rr.stats.component == "R")


def test_rotate_ne_trace_to_rt_theta_quarter_turn():
    """
    At theta==pi/2, theta_p==0: cos(theta)==0, sin(theta)==1, cos(theta_p)==1,
    sin(theta_p)==0, so Eq. 9 reduces to TT==-EN and Eq. 10 reduces to RR==NE
    """
    tr_ee, tr_ne = _trace([1., 2., 3.]), _trace([10., 20., 30.])
    tr_en, tr_nn = _trace([100., 200., 300.]), _trace([4., 5., 6.])

    tr_rr, tr_tt = rotate_ne_trace_to_rt(tr_ee=tr_ee, tr_ne=tr_ne, tr_en=tr_en,
                                         tr_nn=tr_nn, theta=np.pi / 2,
                                         theta_p=0.)

    assert(tr_tt.data == pytest.approx(-tr_en.data))
    assert(tr_rr.data == pytest.approx(tr_ne.data))


def test_rotate_rt_adjsrc_to_ne_invalid_component():
    """Only 'R' or 'T' component adjoint sources can be rotated"""
    tr = _trace([1., 2., 3.])
    tr.stats.component = "Z"
    with pytest.raises(AssertionError):
        rotate_rt_adjsrc_to_ne(tr=tr, theta=0., theta_p=0.)


def test_rotate_rt_adjsrc_to_ne_t_component():
    """
    At theta==theta_p==0, Eq. 16 reduces to: EE==input, EN==NE==NN==0
    """
    tr = _trace([1., 2., 3.])
    tr.stats.component = "T"

    tr_ee, tr_ne, tr_en, tr_nn = rotate_rt_adjsrc_to_ne(tr=tr, theta=0.,
                                                        theta_p=0.)

    assert(tr_ee.data == pytest.approx(tr.data))
    assert(tr_ne.data == pytest.approx(np.zeros(3)))
    assert(tr_en.data == pytest.approx(np.zeros(3)))
    assert(tr_nn.data == pytest.approx(np.zeros(3)))


def test_atomic_write_lines_round_trip(tmpdir):
    """Basic sanity check: writes exactly the given lines, nothing more"""
    fid = os.path.join(tmpdir, "state.txt")
    atomic_write_lines(fid, ["a: 1\n", "b: 2\n"])
    with open(fid) as f:
        assert(f.readlines() == ["a: 1\n", "b: 2\n"])


def test_atomic_write_lines_overwrites_existing(tmpdir):
    """A second write should fully replace the first, not append/merge"""
    fid = os.path.join(tmpdir, "state.txt")
    atomic_write_lines(fid, ["a: 1\n"])
    atomic_write_lines(fid, ["b: 2\n"])
    with open(fid) as f:
        assert(f.readlines() == ["b: 2\n"])


def test_atomic_write_lines_no_leftover_tmp_file(tmpdir):
    """
    The temp file used internally should never remain after a successful
    write -- only the real target file should be left in the directory
    """
    fid = os.path.join(tmpdir, "state.txt")
    atomic_write_lines(fid, ["a: 1\n"])
    assert(glob(os.path.join(tmpdir, "*")) == [fid])


def test_atomic_write_lines_preserves_original_on_failure(tmpdir, monkeypatch):
    """
    If writing fails partway through, the ORIGINAL file must be left
    untouched (not truncated/corrupted) -- this is the whole point of
    writing to a temp file first rather than overwriting in place
    """
    fid = os.path.join(tmpdir, "state.txt")
    atomic_write_lines(fid, ["original: 1\n"])

    def _boom(*args, **kwargs):
        raise OSError("simulated crash mid-write")

    monkeypatch.setattr(os, "replace", _boom)
    with pytest.raises(OSError):
        atomic_write_lines(fid, ["new: 2\n"])

    # Original file is untouched, and no stray temp file left behind
    with open(fid) as f:
        assert(f.readlines() == ["original: 1\n"])
    assert(glob(os.path.join(tmpdir, "*")) == [fid])


def test_rotate_rt_adjsrc_to_ne_r_component():
    """
    At theta==theta_p==0, Eq. 18 reduces to: NN==input, EE==NE==EN==0
    """
    tr = _trace([1., 2., 3.])
    tr.stats.component = "R"

    tr_ee, tr_ne, tr_en, tr_nn = rotate_rt_adjsrc_to_ne(tr=tr, theta=0.,
                                                        theta_p=0.)

    assert(tr_nn.data == pytest.approx(tr.data))
    assert(tr_ee.data == pytest.approx(np.zeros(3)))
    assert(tr_ne.data == pytest.approx(np.zeros(3)))
    assert(tr_en.data == pytest.approx(np.zeros(3)))


@pytest.mark.parametrize("seconds,expected", [
    (0, "0s"), (3.4, "3s"), (59.6, "1m00s"), (125, "2m05s"),
    (3600, "1h00m"), (3720, "1h02m"), (90000, "25h00m"),
])
def test_msg_duration(seconds, expected):
    """Durations in the log are as short as their scale allows"""
    from seisflows.tools import msg
    assert(msg.duration(seconds) == expected)
