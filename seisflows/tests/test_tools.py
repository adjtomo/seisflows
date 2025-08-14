"""
Test any of the utility functions defined in the Tools directory

.. note::

    The Model files in ./tests/
"""
import os
import pytest
import numpy as np
from seisflows import ROOT_DIR
from seisflows.tools.specfem_model import Model
from seisflows.tools.config import custom_import
from seisflows.tools.specfem import read_fortran_binary


TEST_DIR = os.path.join(ROOT_DIR, "tests")
TEST_MODEL = os.path.join(TEST_DIR, "test_data", "test_tools", 
                          "test_specfem_model")


@pytest.fixture
def test_model_serial():
    """
    Provide a Model object to be manipulated in serial

    :rtype: seisflows.specfem_model.Model
    :return: Model object
    """
    m = Model(path=TEST_MODEL, parameters=["c11", "c22", "c33"], parallel=False)
    return m


@pytest.fixture
def test_model_parallel():
    """
    Provide a Model object to be manipulated in prallel

    :rtype: seisflows.specfem_model.Model
    :return: Model object
    """
    m = Model(path=TEST_MODEL, parameters=["c11", "c22", "c33"], parallel=True)
    return m

@pytest.fixture
def test_model_other_serial(tmpdir, test_model_serial):
    """
    Provide a different Model that can be applied with original model

    :rtype: seisflows.specfem_model.Model
    :return: Model object
    """
    test_model_serial.apply(actions=["*"], values=[2], export_to=tmpdir)
    m = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=False)

    return m


@pytest.fixture
def test_model_other_parallel(tmpdir, test_model_parallel):
    """
    Provide a different Model that can be applied with original model

    :rtype: seisflows.specfem_model.Model
    :return: Model object
    """
    test_model_parallel.apply(actions=["*"], values=[2], 
                              export_to=tmpdir)
    m = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=True)

    return m

def test_model_read(test_model_serial):
    """
    Check that model values are read in correctly and accessible in a way we 
    expect
    """
    parameters = ["c11", "c22", "c33"]
    assert(test_model_serial.fmt == ".bin")
    assert(test_model_serial.nproc == 3)
    assert(len(test_model_serial.filenames) == \
           test_model_serial.nproc * len(parameters))


def test_model_apply_wo_other_serial(tmpdir, test_model_serial):
    """
    Test out the main apply function without involving another model
    """
    m = test_model_serial
    arr = m.read(m.filenames[0])
    
    # Make sure there are values in the model
    assert(arr.max() != 0.)

    # Multiply the input array by 0 and export, check that this works as expect
    m.apply(actions=["*"], values=[0], export_to=tmpdir)

    # Read the output model and check that the values are all 0
    m_out = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=False)
    arr_out = m_out.read(m_out.filenames[0])
    assert((arr_out == 0).all())


def test_model_apply_chain_wo_other_serial(tmpdir, test_model_serial):
    """
    Test chaining of multiple actions without other model
    """
    m = test_model_serial
    arr = m.read(m.filenames[0])
    assert(arr.max() != 0.)

    # Test out the addition/subtraction feature and chaining actions
    actions = ["*", "+", "+"] 
    values = [0, -1, 3]
    check_output = sum(values) 

    m.apply(actions=actions, values=values, export_to=tmpdir)
    m = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=False)
    arr_out = m.read(m.filenames[0])

    assert((arr_out == check_output).all())


def test_model_apply_wo_other_parallel(tmpdir, test_model_parallel):
    """
    Same as above but for parallel
    """
    m = test_model_parallel
    arr = m.read(m.filenames[0])
    
    # Make sure there are values in the model
    assert(arr.max() != 0.)

    # Multiply the input array by 0 and export, check that this works as expect
    m.apply(actions=["*"], values=[0], export_to=tmpdir)

    # Read the output model and check that the values are all 0
    m_out = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=True)
    arr_out = m_out.read(m_out.filenames[0])
    assert((arr_out == 0).all())


def test_model_apply_chain_wo_other_parallel(tmpdir, test_model_parallel):
    """
    Test chaining of multiple actions without other model
    """
    m = test_model_parallel
    arr = m.read(m.filenames[0])
    assert(arr.max() != 0.)

    # Test out the addition/subtraction feature and chaining actions
    actions = ["*", "+", "+"] 
    values = [0, -1, 3]
    check_output = sum(values) 

    m.apply(actions=actions, values=values, export_to=tmpdir)
    m = Model(path=tmpdir, parameters=["c11", "c22", "c33"], parallel=True)
    arr_out = m.read(m.filenames[0])

    assert((arr_out == check_output).all())


def test_model_apply_w_other_serial(tmpdir, test_model_serial, 
                                    test_model_other_serial):
    """
    Test out the main apply function without involving another model

    TODO: better checks on this test
    """
    # Run through all the actions and check if values changed
    arr = test_model_serial.read(test_model_serial.filenames[0])
    arr_other = test_model_other_serial.read(
        test_model_other_serial.filenames[0])
    assert(arr != arr_other).all()
    
    # To ensure this doesn't interfere with the `other` model which is in tmpdir
    export_to = os.path.join(tmpdir, "exported_model")

    for action in test_model_serial.acceptable_actions:
        print(action)
        test_model_serial.apply(actions=[action], other=test_model_other_serial,
                                export_to=export_to)
        m = Model(path=export_to, parameters=["c11", "c22", "c33"], 
                  parallel=False)
        arr_check = m.read(m.filenames[0])
        assert(arr[0] != arr_check[0])


def test_model_apply_w_other_parallel(tmpdir, test_model_parallel, 
                                      test_model_other_parallel):
    """
    Test out the main apply function without involving another model

    TODO: better checks on this test
    """
    # Run through all the actions and check if values changed
    arr = test_model_parallel.read(test_model_parallel.filenames[0])
    arr_other = test_model_other_parallel.read(
        test_model_other_parallel.filenames[0])
    assert(arr != arr_other).all()
    
    # To ensure this doesn't interfere with the `other` model which is in tmpdir
    export_to = os.path.join(tmpdir, "exported_model")

    for action in test_model_parallel.acceptable_actions:
        print(action)
        test_model_parallel.apply(actions=[action], 
                                  other=test_model_other_parallel,
                                  export_to=export_to)
        m = Model(path=export_to, parameters=["c11", "c22", "c33"], 
                  parallel=False)
        arr_check = m.read(m.filenames[0])
        assert(arr[0] != arr_check[0])

def test_model_dot_serial(test_model_serial, test_model_other_serial):
    """
    Test out the main apply function without involving another model
    """
    # Brute force compute the dot product manually so we know what the right
    # answer is to compare to
    model_vector, other_model_vector = np.array([]), np.array([])
    for fid in test_model_serial.filenames:
        model_vector = np.append(model_vector, read_fortran_binary(fid))
    for fid in test_model_other_serial.filenames:
        other_model_vector = np.append(other_model_vector, 
                                       read_fortran_binary(fid))
    dot_product_to_check = np.dot(model_vector, other_model_vector)

    # Generate dot product with Model
    dot_product = test_model_serial.dot(test_model_other_serial)

    # There will be some floating point rounding issues but they should be 
    # more or less the same
    assert(dot_product == pytest.approx(dot_product_to_check, 1E-4))

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
