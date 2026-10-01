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
from seisflows.preprocess.default import (Default, read, write,
                                          initialize_adjoint_traces)
from seisflows.preprocess.pyaflowa import Pyaflowa


TEST_DATA = os.path.join(ROOT_DIR, "tests", "test_data", "test_preprocess")
TEST_SOLVER = os.path.join(ROOT_DIR, "tests", "test_data", "test_solver")


def test_default_read():
    """
    Test that we can read observed/synthetic traces with the preprocess mod.

    .. note::
        `Default.read`/`.write`/`.initialize_adjoint_traces` were pulled out
        of the class and now live as module-level functions in
        `seisflows.preprocess.default` (still used internally by `Default`
        and by other preprocess-adjacent workflows, e.g.
        `workflow.noise_inversion`)
    """
    # If new data formats are added to preprocess, they need to be tested
    tested_data_formats = ["ASCII", "SU", "SAC"]

    preprocess = Default()
    assert(set(tested_data_formats) ==
           set(preprocess._obs_acceptable_data_formats))

    st1 = read(os.path.join(TEST_DATA, "AA.S0001.BXY.semd"),
              data_format="ascii")

    st2 = read(os.path.join(TEST_DATA, "Uy_file_single_d.su"),
              data_format="su")

    st3 = read(os.path.join(TEST_DATA, "AA.S0001.BXY.sac"),
              data_format="sac")

    assert(st1[0].stats.npts == st2[0].stats.npts)
    assert(st3[0].stats.npts == st2[0].stats.npts)


def test_default_write(tmpdir):
    """
    Test that we can write synthetic waveforms to formats that SPECFEM recognizes
    """
    # If new data formats supported by SPECFEM are added to preprocess,
    # they need to be tested
    tested_data_formats = ["ASCII", "SU"]

    preprocess = Default()
    assert(set(tested_data_formats) ==
           set(preprocess._syn_acceptable_data_formats))

    st1 = read(os.path.join(TEST_DATA, "AA.S0001.BXY.semd"),
              data_format="ascii")

    write(st1, fid=os.path.join(tmpdir, "test_stream_ascii"),
         data_format="ASCII")
    write(st1, fid=os.path.join(tmpdir, "test_stream_su"), data_format="SU")

    assert(os.path.exists(os.path.join(tmpdir, "test_stream_ascii")))
    assert(os.path.exists(os.path.join(tmpdir, "test_stream_su")))


def test_default_initialize_adjoint_traces(tmpdir):
    """
    Make sure we can write empty adjoint sources expected by SPECFEM
    """
    data_filenames = glob(os.path.join(TEST_DATA, "*semd"))
    initialize_adjoint_traces(data_filenames=data_filenames, fmt="ASCII",
                              path_out=tmpdir)

    data_filenames = glob(os.path.join(TEST_DATA, "*su"))
    initialize_adjoint_traces(data_filenames=data_filenames, fmt="SU",
                              path_out=tmpdir)

    assert(len(glob(os.path.join(tmpdir, "*"))) == 2)
    for fid in glob(os.path.join(tmpdir, "*")):
        assert(fid.endswith(".adj"))


def test_default_quantify_misfit(tmpdir):
    """
    Quantify misfit with some example data
    """
    preprocess = Default(syn_data_format="ascii", obs_data_format="ascii",
                         unit_output="disp", misfit="waveform",
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
    # preprocess.syn_data_format = "SU"
    # preprocess.obs_data_format = "SU"
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

    .. note::
        Discovering `_station_codes`/`_source_names` ahead of time used to be
        `setup()`'s job, but that responsibility has moved to the `Solver`
        module (the workflow now passes an explicit `source_name` into
        `quantify_misfit()` per call; see `workflow.forward.
        evaluate_objective_function`). `setup()` today only builds the
        shared Pyatoa `Config` template and the scratch directory structure.
    """
    pyaflowa = Pyaflowa(
        workdir=tmpdir,
        path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
        path_solver=os.path.join(TEST_SOLVER, "mainsolver"),
        source_prefix="SOURCE",
        ntask=2,
        components="Y",
    )
    assert(pyaflowa._config is None)

    pyaflowa.setup()

    assert(pyaflowa._config is not None)
    assert(pyaflowa._config.component_list == ["Y"])
    for pathname in ["scratch", "_logs", "_tmplogs", "_datasets", "_figures",
                     "_preproc_output"]:
        assert(os.path.exists(pyaflowa.path[pathname]))


def test_pyaflowa_setup_config(tmpdir):
    """
    Test that per-evaluation Config objects (used to control
    `quantify_misfit`) are correctly derived from the shared Config template
    built by `setup()`
    """
    pyaflowa = Pyaflowa(
        workdir=tmpdir,
        path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
        path_solver=TEST_SOLVER, source_prefix="SOURCE", ntask=1,
        data_case="synthetic", components="Y", fix_windows="ITER",
    )
    pyaflowa.setup()

    # This mirrors what `quantify_misfit()` does internally to derive a
    # per-source, per-evaluation Config from the shared template
    config = pyaflowa._config.copy()
    config.event_id = "001"
    config.iteration = 1
    config.step_count = 1
    assert(pyaflowa.ftag(config) == "001_i01s01")


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


@pytest.mark.slow
def test_pyaflowa_line_search(tmpdir):
    """
    Test that the Pyaflowa preprocess class can quantify misfit over the course
    of a few evaluations (a line search) and run its finalization task
    Essentially an integration test testing the entire preprocessing module
    works as a whole

    .. note::
        Marked `slow`: each `quantify_misfit()` call here runs real
        Pyflex/Pyadjoint window-picking and misfit measurement (not mocked),
        so this test alone takes several minutes.
    """
    pyaflowa = Pyaflowa(
        workdir=tmpdir,
        path_specfem_data=os.path.join(TEST_SOLVER, "mainsolver", "DATA"),
        path_output=os.path.join(tmpdir, "output"),
        path_solver=TEST_SOLVER, source_prefix="SOURCE", ntask=2,
        data_case="synthetic", obs_data_format="ASCII", components="Y",
        fix_windows="ITER",
        export_datasets=True, export_figures=True, export_log_files=True,
    )
    pyaflowa.setup()
    unix.mkdir(pyaflowa.path.output)  # usually done by other modules setup
    # `setup()` no longer discovers source names ahead of time (see
    # `test_pyaflowa_setup`); the calling workflow already knows them via
    # `solver.source_names`, which for `ntask=2` with this SOURCE-prefixed
    # test fixture are "001" and "002" (see `test_solver.py::test_source_names`)
    save_residuals = os.path.join(tmpdir, "residuals.txt")
    for source_name in ["001", "002"]:
        for step_count in range(3):
            # Ignore any outputs, just want to run misfit quantification
            # misfit will not be reducing but thats okay
            pyaflowa.quantify_misfit(source_name=source_name,
                                     iteration=1, step_count=step_count,
                                     save_residuals=save_residuals,
                                     save_adjsrcs=tmpdir)

    pyaflowa.finalize()

    # `quantify_misfit()` overwrites (not appends to) `save_residuals` on
    # every call, so this is just the last evaluation's normalized misfit
    final_misfit = float(open(save_residuals).read().strip())
    assert(final_misfit == pytest.approx(0.161, rel=1e-2))

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
