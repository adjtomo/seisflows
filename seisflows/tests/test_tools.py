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
    assert(m.fmt == ".bin")
    assert(m.nproc == 48)
    assert(len(m.filenames) == m.nproc * len(parameters))


def test_model_apply_wo_other(tmpdir):
    """
    Test out the main apply function without involving another model
    """
    m = Model(path=TEST_MODEL, parameters=["c11", "c22", "c33"], parallel=False)
    arr = m.read(m.filenames[0])
    
    # Make sure there are values in the model
    assert(arr.max() != 0.)

    # Multiply the input array by 0 and export, check that this works as expect
    m.apply(actions=["multiply"], values=[0], export_to=tmpdir)

    # Read the output model and check that the values are all 0
    m_out = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=False)
    arr_out = m_out.read(m_out.filenames[0])
    assert((arr_out == 0).all())


def test_model_apply_chain_wo_other(tmpdir):
    """
    Test chaining of multiple actions without other model
    """
    m = Model(path=TEST_MODEL, parameters=["c11", "c22", "c33"], parallel=False)
    arr = m.read(m.filenames[0])
    assert(arr.max() != 0.)

    # Test out the addition/subtraction feature and chaining actions
    actions = ["multiply", "add", "add"] 
    values = [0, -1, 3]
    check_output = sum(values) 

    m_out.apply(actions=actions, values=values, export_to=tmpdir)
    m_out = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=False)
    arr_out = m_out.read(m_out.filenames[0])

    assert((arr_out == check_output).all())


def test_model_loop_parallel(tmpdir):
    """
    Same as single function but process in parallel with concurrent futures
    """
    m = Model(path=TEST_MODEL, parameters=["c11", "c22", "c33"], parallel=True)
    arr = m.read(m.filenames[0])
    # Make sure there are values in the model
    assert(arr.max() != 0.)

    # Multiply the input array by 0 and export, check that this works as expect
    m.apply(action="multiply", value=0, export_to=tmpdir)

    # Read the output model and check that the values are all 0
    m_out = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=True)
    arr_out = m_out.read(m_out.filenames[0])
    assert((arr_out == 0).all())

    # Test out the addition feature 
    m_out.apply(action="add", value=-1, export_to=tmpdir)
    m_out = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=True)
    arr_out = m_out.read(m_out.filenames[0])
    assert((arr_out == -1.).all())

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
