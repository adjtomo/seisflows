"""
Test the ability of the Solver module to interact with various versions of
SPECFEM
"""
import os
import pytest
from glob import glob
from seisflows import ROOT_DIR
from seisflows.tools.config import set_task_id
from seisflows.solver.specfem import Specfem


TEST_DATA = os.path.join(ROOT_DIR, "tests", "test_data", "test_solver")


def test_source_names():
    """
    Check that source names are established correctly
    """
    # Establish solver that looks for CMTSOLUTIONS with NTASK==1
    sources = os.path.join(TEST_DATA, "sources")
    solver = Specfem(path_specfem_data=sources, source_prefix="CMTSOLUTION")
    assert(len(solver.source_names) == 1)

    # Set NTASK==2 to grab both source files
    solver = Specfem(path_specfem_data=sources, source_prefix="CMTSOLUTION",
                     ntask=2)
    source_names = glob(os.path.join(sources, "CMTSOLUTION*"))
    source_names = [_.split("_")[-1] for _ in source_names]

    assert(source_names == solver.source_names)


def test_initialize_working_directory(tmpdir):
    """
    Test that data filenames are returned correctly
    """
    specfem_data = os.path.join(TEST_DATA, "mainsolver", "DATA")
    specfem_bin = os.path.join(TEST_DATA, "mainsolver", "bin")

    solver = Specfem(path_specfem_data=specfem_data,
                     path_specfem_bin=specfem_bin,
                     source_prefix="SOURCE", workdir=tmpdir
                     )
    assert(not os.path.exists(solver.path.mainsolver))

    # Set the environment task id so that the logger doesn't throw warnings
    # about not finding the task id
    set_task_id(0)

    # Generate the required directory structure
    solver._initialize_working_directory(
        cwd=os.path.join(solver.path.scratch, "001")
    )

    # Simple checks to make sure the directory structure was set up properly
    assert(os.path.islink(solver.path.mainsolver))
    assert(os.path.exists(solver.cwd))
    assert(glob(os.path.join(solver.cwd, "*")))
    event_fid = os.path.join(solver.cwd, "DATA", "SOURCE")
    assert(os.path.islink(event_fid))
    event_line = open(event_fid).readlines()[0].strip()
    assert(event_line == "## Source 1")


def test_run_binary(tmpdir):
    """
    Just run a known and intentially incorrect binary with the run binary
    function to check that we can, and that error catching is working
    """
    solver = Specfem()
    solver._run_binary(executable="echo hello world",
                       stdout=os.path.join(tmpdir, "log.txt"))

    # Executables that don't exist will not run
    with pytest.raises(SystemExit):
        solver._run_binary(executable="gobbledigook",
                           stdout=os.path.join(tmpdir, "log.txt"))

    # Executables that do exist but error out (e.g., with invalid options) will
    # also throw an error
    with pytest.raises(SystemExit):
        solver._run_binary(executable="ls -//daflkjeaf",
                           stdout=os.path.join(tmpdir, "log.txt"))
