"""
Test the ability of the Solver module to interact with various versions of
SPECFEM
"""
import os
import pytest
import numpy as np
from glob import glob
from pyasdf import ASDFDataSet
from seisflows import ROOT_DIR
from seisflows.tools import unix
from seisflows.preprocess.default import Default
from seisflows.preprocess.pyaflowa import Pyaflowa


TEST_DATA = os.path.join(ROOT_DIR, "tests", "test_data", "test_preprocess")
TEST_SOLVER = os.path.join(ROOT_DIR, "tests", "test_data", "test_solver")


def test_read():
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


def test_write(tmpdir):
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


def test_initialize_adjoint_traces(tmpdir):
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


def test_quantify_misfit(tmpdir):
    """
    Quantify misfit with some example data
    """
    preprocess = Default(data_format="ascii", misfit="waveform",
                         adjoint="waveform", path_preprocess=tmpdir)
    preprocess.setup()

    data_filenames = glob(os.path.join(TEST_DATA, "*semd"))
    preprocess.quantify_misfit(
        observed=data_filenames, synthetic=data_filenames,
        save_residuals=os.path.join(tmpdir, "residuals_ascii"),
        save_adjsrcs=tmpdir
    )

    preprocess.data_format = "SU"
    data_filenames = glob(os.path.join(TEST_DATA, "*su"))
    preprocess.quantify_misfit(
        observed=data_filenames, synthetic=data_filenames,
        save_residuals=os.path.join(tmpdir, "residuals_su"),
        save_adjsrcs=tmpdir
    )

    assert(len(glob(os.path.join(tmpdir, "*"))) == 4)
    residuals = open(os.path.join(tmpdir, "residuals_ascii")).readlines()
    assert(len(residuals) == 1)
    assert(float(residuals[0]) == 0)


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

    assert(len(pyaflowa._station_codes) == 5)
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


def test_pyaflowa_quantify_misfit_station(tmpdir):
    """
    Check that the function to quantify misfit that should be run in parallel
    works as a serial job
    """
    pyaflowa = Pyaflowa(
        workdir=tmpdir,
        path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
        path_solver=TEST_SOLVER, source_prefix="SOURCE", ntask=2,
        data_case="synthetic", components="Y",
    )
    pyaflowa.setup()
    config = pyaflowa._setup_quantify_misfit(source_name="001", iteration=1,
                                             step_count=1)
    misfit, nwin = pyaflowa._quantify_misfit_station(
        config=config, station_code=pyaflowa._station_codes[0],
        save_adjsrcs=False
    )
    assert(misfit == 33.5304)
    assert(nwin == 8.)


def test_pyaflowa_quantify_misfit_single(tmpdir):
    """
    Test misfit quantification for Pyatoa during a single misfit evaluation.
    Waveform data and source and receiver metadata is exposed from the test data
    directory. Data and synthetics are the same so residuals will be 0. Want
    to check that we can process in parallel and that Pyatoa outputs figures,
    and data.
    """
    pyaflowa = Pyaflowa(
        workdir=tmpdir,
        path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
        path_solver=TEST_SOLVER, source_prefix="SOURCE", ntask=2,
        data_case="synthetic", components="Y",
    )
    pyaflowa.setup()
    for source_name in pyaflowa._source_names:
        save_residuals = os.path.join(tmpdir, f"residuals_{source_name}.txt")
        pyaflowa.quantify_misfit(source_name=source_name,
                                 save_residuals=save_residuals,
                                 save_adjsrcs=tmpdir)

    residuals = np.loadtxt(save_residuals)  # just check one of the file
    assert(residuals == 0.919)

    # Check that windows and adjoint sources were saved to dataset
    nwin = {"001": 45, "002": 48}
    for source_name in pyaflowa._source_names:
        with ASDFDataSet(os.path.join(pyaflowa.path._datasets,
                                      f"{source_name}.h5")) as ds:
            # Pyatoa selects N number windows for each source
            assert(len(ds.auxiliary_data.MisfitWindows.i01.s00.list()) ==
                   nwin[source_name])
            assert(len(ds.auxiliary_data.AdjointSources.i01.s00.list()) == 5)

    # Check that adjoint sources are all zero
    adjsrcs = glob(os.path.join(tmpdir, "*.adj"))
    for adjsrc in adjsrcs:
        data = np.loadtxt(adjsrc)
        assert(data[:, 1].any())  # assert that adjoint sourcse are not zero


def test_pyaflowa_check_fixed_windows(tmpdir):
    """
    Test that misfit window bool returner always returns how we want it to.
    """
    pf = Pyaflowa(fix_windows=True)
    assert(pf._check_fixed_windows(iteration=99, step_count=99)[0])
    pf = Pyaflowa(fix_windows="ITER")
    assert(not pf._check_fixed_windows(iteration=1, step_count=0)[0])
    assert(pf._check_fixed_windows(iteration=1, step_count=1)[0])
    pf = Pyaflowa(fix_windows="ONCE", start=5)
    assert(not pf._check_fixed_windows(iteration=5, step_count=0)[0])
    assert(pf._check_fixed_windows(iteration=5, step_count=1)[0])
    assert(pf._check_fixed_windows(iteration=6, step_count=0)[0])


def test_pyaflowa_finalize(tmpdir):
    """
    Test teardown procedures for the Pyaflowa preprocessing module which
    includes creating an Inspector, condensing PDF files, and exporting
    files to disk.
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
    for source_name in pyaflowa._source_names:
        for step_count in range(3):
            # Ignore any outputs, just want to run misfit quantification
            # misfit will not be reducing but thats okay
            pyaflowa.quantify_misfit(source_name=source_name,
                                     iteration=1,
                                     step_count=step_count)

    pyaflowa.finalize()
    # Just check file count to see that finalize did what it's supposed to do
    # since finalize just moves and collects files
    assert(len(glob(os.path.join(pyaflowa.path.output, "figures", "*"))) == 1)
    assert(len(glob(os.path.join(pyaflowa.path.output, "logs", "*"))) == 6)
    assert(len(glob(os.path.join(pyaflowa.path.output,
                                 "datasets", "*.csv"))) == 2)


# def test_pyaflowa_quantify_misfit_inversion(tmpdir):
#     """
#     Test misfit quantification for Pyatoa but simulating multiple back-to-back
#     evaluations as one would encounter during an inversion This would involve
#     re-using misfit windows throughout the evaluation, and reading in already
#     gathered data from an ASDFDataSet
#     """
#     pyaflowa = Pyaflowa(
#         workdir=tmpdir,
#         path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
#         path_solver=TEST_SOLVER, source_prefix="SOURCE", ntask=1,
#         data_case="synthetic", components="Y", fix_windows="ITER",
#     )
#     pyaflowa.setup()
#     source_name = pyaflowa._source_names[0]
#     for step_count in range(2):
#         # Ignore any outputs, just want to run misfit quantification
#         # misfit will not be reducing but thats okay
#         pyaflowa.quantify_misfit(source_name=source_name,
#                                  iteration=1,
#                                  step_count=step_count)
#
#     # Check that correct number of PDFs have been made
#     assert(len(glob(os.path.join(pyaflowa.path._figures, "*pdf"))) == 4)
#
#     # Check datasets for correct formatting of auxiliary data and rand vals
#     fid = os.path.join(pyaflowa.path._datasets, f"{source_name}.h5")
#     with ASDFDataSet(fid, mode="r") as ds:
#         assert(len(ds.waveforms.list()) == 5)
#         sta_0 = ds.waveforms[ds.waveforms.list()[0]]
#         assert(len(sta_0.list()) == 5)  # 1 observed, 4 synthetics
#         assert(len(ds.auxiliary_data.AdjointSources.i01) == 2)
#         assert(len(ds.auxiliary_data.MisfitWindows.i01) == 2)
#         assert(len(ds.auxiliary_data.MisfitWindows.i01.s01.list()) == 45)
#         adjsrc = ds.auxiliary_data.AdjointSources.i01.s00.AA_S000004_BXY
#         misfit = adjsrc.parameters["misfit"]
#         assert(misfit == pytest.approx(18.3167, 3))
