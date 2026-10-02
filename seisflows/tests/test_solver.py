"""
Test the ability of the Solver module to interact with various versions of
SPECFEM
"""
import os
import shutil
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
    assert(not os.path.exists(solver.path._mainsolver))

    # Set the environment task id so that the logger doesn't throw warnings
    # about not finding the task id
    set_task_id(0)

    # Generate the required directory structure
    solver._initialize_working_directory(
        cwd=os.path.join(solver.path.scratch, "001")
    )

    # Simple checks to make sure the directory structure was set up properly
    assert(os.path.islink(solver.path._mainsolver))
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


def _write_model_files(path, parameters, nproc=2):
    """Create empty model files named the way SPECFEM names them"""
    for par in parameters:
        for i in range(nproc):
            open(os.path.join(path, f"proc{i:06d}_{par}.bin"), "w").close()


def test_check_model_files(tmpdir):
    """A complete set of model files passes the cheap structural check"""
    solver = Specfem(materials="ELASTIC")
    _write_model_files(str(tmpdir), solver._parameters)
    solver.check_model_files(path=str(tmpdir))  # does not exit


def test_check_model_files_does_not_read_values(tmpdir):
    """
    The point of the structural check is that the master job can run it: the
    files here are empty, so reading any model values out of them would fail
    """
    solver = Specfem()
    _write_model_files(str(tmpdir), solver._parameters)
    solver.check_model_files(path=str(tmpdir))


@pytest.mark.parametrize("problem", ["no_path", "missing_parameter",
                                     "incomplete_write"])
def test_check_model_files_catches_structural_problems(tmpdir, problem):
    solver = Specfem(materials="ELASTIC")  # parameters vp, vs
    path = str(tmpdir)

    if problem == "no_path":
        path = os.path.join(path, "does_not_exist")
    elif problem == "missing_parameter":
        _write_model_files(path, ["vp"])  # no vs
    elif problem == "incomplete_write":
        _write_model_files(path, ["vp"], nproc=2)
        _write_model_files(path, ["vs"], nproc=1)  # fewer files than vp

    with pytest.raises(SystemExit):
        solver.check_model_files(path=path)


def test_validate_model_values_passes_for_real_model(tmpdir):
    solver = Specfem(materials="ELASTIC")
    # The test data directory holds the model in several formats at once,
    # which `Model` can't auto-detect; a real model directory has one
    src = os.path.join(ROOT_DIR, "tests", "test_data", "test_tools",
                       "test_file_formats")
    for fid in glob(os.path.join(src, "*.bin")):
        shutil.copy(fid, str(tmpdir))
    solver.validate_model_values(path=str(tmpdir))  # does not raise


def test_validate_model_values_raises_instead_of_exiting(tmpdir):
    """
    Unlike `check_model_values` this raises, so that the workflow can record
    why a compute task failed instead of the task just exiting
    """
    solver = Specfem()
    with pytest.raises(AssertionError):
        solver.validate_model_values(path=os.path.join(str(tmpdir), "nope"))
