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
from seisflows.tools.state import State


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
    m_new = Model(path=os.path.join(tmpdir, "test.npz"))
    assert(m_new.ngll[0] == m.ngll[0])
    assert(m_new.fmt == m.fmt)

    # Check that writing fortran binary works
    m.write(path=tmpdir)
    assert(len(glob(os.path.join(tmpdir, f"*{m.fmt}"))) == 2)

    # Check that we can instantiate a model from an input vector
    m = Model()
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
    

def test_state_arr_2_str(tmpdir):
    """
    Test the state checkpointing system conversion system for arrays to strings
    with a few test cases
    """
    state = State(fid=os.path.join(tmpdir, "sfstate.json"))

    # Check arr2str and str2arr functionality
    arr = [0,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,0,0,1,1,0]
    str_ = "0,4-6,9,14-17,20"
    assert(state._arr_to_str(arr) == str_)
    assert(state._str_to_arr(str_, ntasks=len(arr)) == arr)    

    arr = [0] * 100
    str_ = "0-99"
    assert(state._arr_to_str(arr) == str_)
    assert(state._str_to_arr(str_, ntasks=len(arr)) == arr)    

    arr = [1] * 100
    str_ = ""
    assert(state._arr_to_str(arr) == str_)
    assert(state._str_to_arr(str_, ntasks=len(arr)) == arr)    

    arr = [0]
    str_  = "0"
    assert(state._arr_to_str(arr) == str_)
    assert(state._str_to_arr(str_, ntasks=len(arr)) == arr)    

    arr = [1]
    str_ = ""
    assert(state._arr_to_str(arr) == str_)
    assert(state._str_to_arr(str_, ntasks=len(arr)) == arr)    


def test_state_read_write_single(tmpdir):
    """
    Test the state checkpointing system for reading and writing to disk
    """
    state = State(fid=os.path.join(tmpdir, "sfstate.json"))
    state("test_function", 10)
    state("test_function_2", 33)
    state("test_function_3", 1000)
    assert(state("test_function", 10) == "0-9")
    assert(state("test_function_2", 33) == "0-32")
    assert(state("test_function_3", 1000) == "0-999")

def test_state_done_single(tmpdir):
    """
    Test the Done feature of completing tasks one by one
    """
    state = State(fid=os.path.join(tmpdir, "sfstate.json"))
    state("test_function", 10)

    for val in [1, 3, 6, 7, 9]:
        state.done("test_function", val)
    assert(state("test_function", 10) == "0,2,4-5,8")

def test_state_done_parallel(tmpdir):
    """
    Test the lock file capabilities by completing tasks in parallel with
    concurrent futures to mimic a parallel environment
    """
    from concurrent.futures import ProcessPoolExecutor, wait
    ntasks = 500
    state = State(fid=os.path.join(tmpdir, "sfstate.json"))
    state("test_function", ntasks)

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(state.done, "test_function", i) 
                   for i in range(0, ntasks)]
    wait(futures)
    for future in futures:
        future.result()
    
    
    assert(state("test_function", ntasks) == "")
