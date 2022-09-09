"""
Test the ability of the Solver module to interact with various versions of
SPECFEM
"""
import os
import numpy as np
import pytest
from glob import glob
from seisflows import ROOT_DIR
from seisflows.tools import unix
from seisflows.preprocess.default import Default
from seisflows.preprocess.pyaflowa import Pyaflowa


TEST_DATA = os.path.join(ROOT_DIR, "tests", "test_data", "test_preprocess")
TEST_SOLVER = os.path.join(ROOT_DIR, "tests", "test_data", "test_solver")


def test_default_read():
    """
    Test that we can read SPECFEM generated synthetics with the preprocess mod.
    """
    # If new data formats are added to preprocess, they need to be tested
    tested_data_formats = ["ASCII", "SU"]
    preprocess = Default()
    assert(set(tested_data_formats) == set(preprocess._acceptable_data_formats))

    preprocess = Default(data_format="ascii")
    st1 = preprocess.read(os.path.join(TEST_DATA, "AA.S0001.BXY.semd"))

    preprocess = Default(data_format="su")
    st2 = preprocess.read(os.path.join(TEST_DATA, "Uy_file_single_d.su"))

    assert(st1[0].stats.npts == st2[0].stats.npts)


def test_default_write(tmpdir):
    """
    Make sure we can write both data formats
    """
    # If new data formats are added to preprocess, they need to be tested
    tested_data_formats = ["ASCII", "SU"]
    preprocess = Default()
    assert(set(tested_data_formats) == set(preprocess._acceptable_data_formats))

    preprocess = Default(data_format="ascii")
    st1 = preprocess.read(os.path.join(TEST_DATA, "AA.S0001.BXY.semd"))
    preprocess.write(st1, fid=os.path.join(tmpdir, "test_stream_ascii"))

    preprocess.data_format = "SU"
    preprocess.write(st1, fid=os.path.join(tmpdir, "test_stream_su"))


def test_default_initialize_adjoint_traces(tmpdir):
    """
    Make sure we can write empty adjoint sources expected by SPECFEM
    """
    preprocess = Default(data_format="ascii")
    data_filenames = glob(os.path.join(TEST_DATA, "*semd"))
    preprocess.initialize_adjoint_traces(data_filenames=data_filenames,
                                         output=tmpdir)

    preprocess.data_format = "SU"
    data_filenames = glob(os.path.join(TEST_DATA, "*su"))
    preprocess.initialize_adjoint_traces(data_filenames=data_filenames,
                                         output=tmpdir)

    assert(len(glob(os.path.join(tmpdir, "*"))) == 2)
    for fid in glob(os.path.join(tmpdir, "*")):
        assert(fid.endswith(".adj"))


def test_default_quantify_misfit(tmpdir):
    """
    Quantify misfit with some example data
    """
    preprocess = Default(data_format="ascii", misfit="waveform",
                         adjoint="waveform", path_preprocess=tmpdir,
                         path_solver=TEST_SOLVER, source_prefix="SOURCE",
                         ntask=2,
                         )
    preprocess.setup()

    preprocess.quantify_misfit(
        source_name="001",
        save_residuals=os.path.join(tmpdir, "residuals_ascii"),
        save_adjsrcs=tmpdir
    )

    # !!! throws a segy error because data are not in the right format
    # preprocess.data_format = "SU"
    # preprocess.quantify_misfit(
    #     source_name="001",
    #     save_residuals=os.path.join(tmpdir, "residuals_su"),
    #     save_adjsrcs=tmpdir
    # )

    assert(len(glob(os.path.join(tmpdir, "*"))) == 3)
    residuals = open(os.path.join(tmpdir, "residuals_ascii")).readlines()
    assert(len(residuals) == 2)
    assert(float(residuals[0]) == pytest.approx(0.0269, 3))


def test_pyaflowa_setup(tmpdir):
    """
    Test setup procedure for SeisFlows which internalizes some workflow
    information that is crucial for later tasks
    """
    pyaflowa = Pyaflowa(
        workdir=tmpdir,
        path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
        path_solver=os.path.join(TEST_SOLVER, "mainsolver"),
        source_prefix="SOURCE",
        ntask=2,
        components="Y",
    )

    assert(pyaflowa._station_codes == [])
    assert(pyaflowa._source_names == [])

    pyaflowa.setup()

    assert(len(pyaflowa._station_codes) == 2)
    assert(pyaflowa._station_codes[0] == "AA.S000000.*.*")
    assert(len(pyaflowa._source_names) == pyaflowa._ntask)
    assert(pyaflowa._source_names[0] == "001")
    assert(pyaflowa._config.component_list == ["Y"])


def test_pyaflowa_setup_quantify_misfit(tmpdir):
    """
    Test Config setup that is used to control `quantify_misfit` function
    """
    pyaflowa = Pyaflowa(
        workdir=tmpdir,
        path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
        path_solver=TEST_SOLVER, source_prefix="SOURCE", ntask=1,
        data_case="synthetic", components="Y", fix_windows="ITER",
    )
    pyaflowa.setup()
    config = pyaflowa._setup_quantify_misfit(source_name="001", iteration=1,
                                             step_count=1)
    assert(config.eval_tag == "i01s01")
    # Data specific time series values calculated by function
    assert(config.start_pad == 48)
    assert(config.end_pad == 299.94)


def test_pyaflowa_check_fixed_windows():
    """
    Test that misfit window bool returner always returns how we want it to.
    """
    pf = Pyaflowa(fix_windows=True)
    assert (pf._check_fixed_windows(iteration=99, step_count=99)[0])
    pf = Pyaflowa(fix_windows="ITER")
    assert (not pf._check_fixed_windows(iteration=1, step_count=0)[0])
    assert (pf._check_fixed_windows(iteration=1, step_count=1)[0])
    pf = Pyaflowa(fix_windows="ONCE", start=5)
    assert (not pf._check_fixed_windows(iteration=5, step_count=0)[0])
    assert (pf._check_fixed_windows(iteration=5, step_count=1)[0])
    assert (pf._check_fixed_windows(iteration=6, step_count=0)[0])


def test_pyaflowa_line_search(tmpdir):
    """
    Test that the Pyaflowa preprocess class can quantify misfit over the course
    of a few evaluations (a line search) and run its finalization task
    Essentially an integration test testing the entire preprocessing module
    works as a whole
    """
    pyaflowa = Pyaflowa(
        workdir=tmpdir,
        path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
        path_output=os.path.join(tmpdir, "output"),
        path_solver=TEST_SOLVER, source_prefix="SOURCE", ntask=2,
        data_case="synthetic", components="Y", fix_windows="ITER",
        export_datasets=True, export_figures=True, export_log_files=True,
    )
    pyaflowa.setup()
    unix.mkdir(pyaflowa.path.output)  # usually done by other modules setup
    save_residuals = os.path.join(tmpdir, f"residuals.txt")
    for source_name in pyaflowa._source_names:
        for step_count in range(3):
            # Ignore any outputs, just want to run misfit quantification
            # misfit will not be reducing but thats okay
            pyaflowa.quantify_misfit(source_name=source_name,
                                     iteration=1, step_count=step_count,
                                     save_residuals=save_residuals,
                                     save_adjsrcs=tmpdir)

    pyaflowa.finalize()

    # Check that final residuals file is the same
    residuals = np.loadtxt(save_residuals)
    assert(pyaflowa.sum_residuals(residuals) == 6.045)

    # Check that atleast one adjoint sources are not zero
    adjsrcs = glob(os.path.join(tmpdir, "*.adj"))
    data = np.loadtxt(adjsrcs[0])
    assert(data[:, 1].any())  # assert that adjoint sourcse are not zero

    # Just check file count to see that finalize did what it's supposed to do
    # since finalize just moves and collects files
    assert(len(glob(os.path.join(pyaflowa.path.output, "pyaflowa", 
                                 "figures", "*"))) == 1)
    assert(len(glob(os.path.join(pyaflowa.path.output, "pyaflowa", 
                                 "logs", "*"))) == 6)
    assert(len(glob(os.path.join(pyaflowa.path.output, "pyaflowa",
                                 "datasets", "*.csv"))) == 2)
