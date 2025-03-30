"""
Test any of the utility functions defined in the Tools directory
"""
import pytest
import time
from seisflows import ROOT_DIR
from seisflows.tools.specfem_model import Model
from seisflows.tools.config import custom_import


# TEST_DIR = os.path.join(ROOT_DIR, "tests")
# TEST_MODEL = os.path.join(TEST_DIR, "test_data", "test_tools", 
                        #   "test_file_formats")
TEST_MODEL = "/Users/chow/Work/software/seisflows/issue_245/DATABASES_MPI"
TEST_MODEL_OTHER = TEST_MODEL


def test_model_read():
    """
    Check that model values are read in correctly and accessible in a way we 
    expect
    """
    parameters = ["c11", "c22", "c33"]
    m = Model(path=TEST_MODEL, parameters=parameters)
    assert(m.flavor == "3D")
    assert(m.fmt == ".bin")
    assert(m.nproc == 48)
    assert(len(m.filenames) == m.nproc * len(parameters))


def test_model_loop_single(tmpdir):
    """
    Test the core merge and split functions of the Model class to ensure
    they give the same results
    """
    m = Model(path=TEST_MODEL, parameters=["c11", "c22", "c33"], parallel=False)
    arr = m.read(m.filenames[0])
    assert(arr.max() != 0.)

    # Multiply the input array by 0 and export, check that this works as expect
    m.apply(action="multiply", value=0, export_to=tmpdir)
    m_out = Model(path=tmpdir, parameters=["c11", "c22", "c33"], fmt=m.fmt)
    arr_out = m_out.read(m_out.filenames[0])
    assert(arr_out.max() == 0.)

def test_model_loop_parallel(tmpdir):
    """
    Same as single function but process in parallel with concurrent futures
    """
    m = Model(path=TEST_MODEL, parameters=["c11", "c22", "c33"], parallel=True)
    arr = m.read(m.filenames[0])
    assert(arr.max() != 0.)

    m.apply(action="multiply", value=0, export_to=tmpdir)
    m_out = Model(path=tmpdir, parameters=["c11", "c22", "c33"], fmt=m.fmt)
    arr_out = m_out.read(m_out.filenames[0])
    assert(arr_out.max() == 0.)

# def test_model_io(tmpdir):
#     """
#     Check that saving and loading npz file works
#     """
#     m = Model(path=TEST_MODEL, fmt=".bin")

#     m.save(path=os.path.join(tmpdir, "test.npz"))
#     m_new = Model(path=os.path.join(tmpdir, "test.npz"))
#     assert(m_new.ngll[0] == m.ngll[0])
#     assert(m_new.fmt == m.fmt)

#     # Check that writing fortran binary works
#     m.write(path=tmpdir)
#     assert(len(glob(os.path.join(tmpdir, f"*{m.fmt}"))) == 2)


# def test_model_from_input_vector():
#     """Check that we can instantiate a model from an input vector"""
#     m = Model(path=None)
#     m.model = Dict(x=[np.array([-1.2, 1.])])
#     assert(m.nproc == 1)
#     assert(m.ngll == [2])
#     assert(m.parameters == ["x"])


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
