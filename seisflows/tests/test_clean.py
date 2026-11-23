"""
Test cleaning the large files out of a working directory, and the command line
`seisflows clean` that does it
"""
import os
import sys
import pytest
from unittest.mock import patch

from seisflows.seisflows import SeisFlows
from seisflows.tools import clean, msg

# Files of a working directory by what should happen to them. Paths are
# relative to the working directory. `KEPT` are files of the workflow that
# have to survive, `LARGE` the ones that are removed by default
LARGE = {
    "forward_arrays": [
        "scratch/solver/001/OUTPUT_FILES/DATABASES_MPI/"
        "proc000000_save_forward_arrays.bin",
        "scratch/solver/001/OUTPUT_FILES/DATABASES_MPI/"
        "proc000000_absorb_field.bin",
        "scratch/solver/002/OUTPUT_FILES/DATABASES_MPI/"
        "proc000000_save_forward_arrays.bin",
        "scratch/solver/001/OUTPUT_FILES/DATABASES_MPI/"
        "proc000001_reg1_absorb_buffer.bin",
        "scratch/solver/001/ZZ_FWD_ARR/proc000000_save_forward_arrays.bin",
    ],
    "vtk": [
        "scratch/solver/001/OUTPUT_FILES/DATABASES_MPI/proc000000_vp.vtk",
    ],
    "databases": [
        "scratch/solver/001/OUTPUT_FILES/DATABASES_MPI/"
        "proc000000_external_mesh.bin",
        "scratch/solver/001/OUTPUT_FILES/DATABASES_MPI/proc000000_vp.bin",
        "scratch/solver/002/OUTPUT_FILES/DATABASES_MPI/proc000000_vp.bin",
        "scratch/solver/002/DATA/proc000000_vs.bin",  # SPECFEM2D models
    ],
    "kernels": [
        "scratch/eval_grad/kernels/001/proc000000_vp_kernel.bin",
        "scratch/eval_grad/kernels/002/proc000000_vp_kernel.bin",
        "scratch/eval_grad/gradient/proc000000_vp_kernel.bin",
        "scratch/eval_grad/misfit_kernel/proc000000_vp_kernel.bin",
    ],
    "models": [
        "scratch/eval_func/model/proc000000_vp.bin",
        "scratch/eval_grad/model/proc000000_vp.bin",
    ],
    "traces": [
        "scratch/solver/001/traces/syn/XX.S0001.BXZ.semd",
        "scratch/solver/001/traces/adj/XX.S0001.BXZ.adj",
        "scratch/solver/001/SEM/XX.S0001.BXZ.adj",
        "scratch/solver/001/OUTPUT_FILES/XX.S0001.BXZ.semd",
    ],
}
NOT_USUAL = {
    "obs": ["scratch/solver/001/traces/obs/XX.S0001.BXZ.semd"],
    "solver_logs": [
        "scratch/solver/001/OUTPUT_FILES/timestamp000001",
        "scratch/solver/001/OUTPUT_FILES/output_list_sources.txt",
        "scratch/solver/001/OUTPUT_FILES/sr.vtk",
        "scratch/solver/001/OUTPUT_FILES/starttimeloop.txt",
    ],
}
KEPT = [
    "parameters.yaml", "sflog.txt", "sfstate.txt", "logs/1_task0.log",
    "output/MODEL_01/proc000000_vp.bin", "output/solver/001/syn/a.semd",
    "scratch/optimize/m_new.npz", "scratch/optimize/f_new.txt",
    "scratch/optimize/LBFGS/S.dat", "scratch/system/funcs.p",
    "scratch/eval_func/residuals/residuals_001_i01s01.txt",
    "scratch/eval_grad/residuals/residuals_001_i01s00.txt",
    "scratch/solver/001/DATA/Par_file", "scratch/solver/001/DATA/STATIONS",
    "scratch/solver/001/OUTPUT_FILES/output_solver.txt",
    "scratch/solver/001/bin/xspecfem3D",
]


def _touch(path, size=10):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(b"x" * size)


@pytest.fixture
def workdir(tmpdir):
    """A working directory with a file for everything above, and links that
    lead outside of it"""
    workdir = str(tmpdir.mkdir("work"))
    fids = [fid for fids in list(LARGE.values()) + list(NOT_USUAL.values())
            for fid in fids] + KEPT
    for fid in fids:
        _touch(os.path.join(workdir, fid))
    with open(os.path.join(workdir, "scratch/solver/001/DATA/Par_file"),
              "w") as f:
        f.write("LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI"
                "\n")

    # The first solver is linked as the 'mainsolver'
    os.symlink(os.path.join(workdir, "scratch/solver/001"),
               os.path.join(workdir, "scratch/solver/mainsolver"))

    # Things that match what is removed but are not in the working directory,
    # reached from inside of it
    outside = str(tmpdir.mkdir("outside"))
    _touch(os.path.join(outside, "proc000000_vp.bin"))
    _touch(os.path.join(outside, "proc000000_save_forward_arrays.bin"))
    _touch(os.path.join(outside, "kernels", "proc000000_vp_kernel.bin"))
    os.symlink(outside, os.path.join(workdir,
                                     "scratch/solver/002/OUTPUT_FILES_LINKED"))
    os.symlink(os.path.join(outside, "proc000000_vp.bin"),
               os.path.join(workdir, "scratch/solver/002/DATA/proc000009_vp.bin"
                            ))
    os.symlink(outside, os.path.join(workdir, "scratch/solver/003"))
    os.symlink(os.path.join(outside, "kernels"),
               os.path.join(workdir, "scratch/eval_grad/kernels/linked"))
    return workdir


def _found(workdir, **kwargs):
    found = clean.find_files(os.path.join(workdir, "scratch"), **kwargs)
    return {name: sorted(os.path.relpath(fid, workdir) for fid in fids)
            for name, fids in found.items()}


def test_files_are_found_by_category(workdir):
    found = _found(workdir)
    assert(set(found) == set(LARGE))  # the ones that are removed as a rule
    for name, fids in LARGE.items():
        assert(found[name] == sorted(fids)), name


def test_files_that_a_workflow_needs_are_not_found(workdir):
    found = _found(workdir, categories=list(clean.CATEGORIES))
    everything = {fid for fids in found.values() for fid in fids}
    assert(not everything & set(KEPT))
    # Not removed unless asked
    assert(found["obs"] == sorted(NOT_USUAL["obs"]))
    assert(found["solver_logs"] == sorted(NOT_USUAL["solver_logs"]))


def test_nothing_outside_of_scratch_is_found(workdir):
    """Links to a directory or a file elsewhere, which a glob would follow"""
    found = _found(workdir, categories=list(clean.CATEGORIES))
    root = os.path.realpath(workdir) + os.sep
    for fids in found.values():
        for fid in fids:
            assert(not os.path.islink(os.path.join(workdir, fid)))
            assert(os.path.realpath(os.path.join(workdir, fid)).startswith(
                root)), fid
    # The first solver is only counted as the directory that it is
    assert(not any("mainsolver" in fid for fids in found.values()
                   for fid in fids))


def test_a_file_has_one_category_whatever_is_chosen(workdir):
    """Leaving out the forward arrays must not remove them as 'databases',
    which they match as well"""
    forward_arrays = set(_found(workdir)["forward_arrays"])
    for categories in [["databases"], ["databases", "vtk"],
                       ["vtk", "databases", "kernels"]]:
        found = _found(workdir, categories=categories)
        assert(not forward_arrays & {fid for fids in found.values()
                                     for fid in fids})
    assert(set(_found(workdir, categories=["databases"])["databases"]) ==
           set(LARGE["databases"]))


def test_select_categories():
    usual = [name for name, (default, _) in clean.CATEGORIES.items()
             if default]
    assert(clean.select_categories() == usual)
    assert(clean.select_categories(include=["obs"]) == usual + ["obs"])
    assert(clean.select_categories(exclude=["forward_arrays", "kernels"]) ==
           [name for name in usual if name not in ["forward_arrays",
                                                   "kernels"]])
    assert("obs" not in clean.select_categories(include=["obs"],
                                                exclude=["obs"]))
    with pytest.raises(ValueError, match="not a category"):
        clean.select_categories(include=["everything"])


def test_the_directory_structure_survives(workdir):
    dirs = sorted(dirpath for dirpath, _, _ in os.walk(workdir))
    found = clean.find_files(os.path.join(workdir, "scratch"))
    nfiles, nbytes, failed = clean.remove_files(
        found, workdir, record=os.path.join(workdir, "scratch/system/rec.txt"))

    assert(sorted(dirpath for dirpath, _, _ in os.walk(workdir)) == dirs)
    assert(nfiles == sum(len(fids) for fids in LARGE.values()))
    assert(nbytes == 10 * nfiles)
    assert(failed == [])
    for fid in sum(LARGE.values(), []):
        assert(not os.path.exists(os.path.join(workdir, fid)))
    for fid in KEPT + sum(NOT_USUAL.values(), []):
        assert(os.path.exists(os.path.join(workdir, fid))), fid
    # And what was outside
    assert(os.path.exists(os.path.join(workdir, "..", "outside",
                                       "proc000000_vp.bin")))
    assert(os.path.exists(os.path.join(workdir, "..", "outside", "kernels",
                                       "proc000000_vp_kernel.bin")))


def test_what_is_removed_is_recorded(workdir):
    record = os.path.join(workdir, "scratch/system/clean_manifest.txt")
    found = clean.find_files(os.path.join(workdir, "scratch"),
                             categories=["models"])
    clean.remove_files(found, workdir, record=record)
    clean.remove_files(
        clean.find_files(os.path.join(workdir, "scratch"),
                         categories=["vtk"]), workdir, record=record)

    with open(record) as f:
        lines = f.read().splitlines()
    assert([line for line in lines if line.startswith("#")] ==
           [line for line in lines if "cleaned" in line])  # added, not replaced
    rows = sorted(line.split("\t") for line in lines
                  if not line.startswith("#"))
    assert(rows == sorted([
        ["models", "10", "scratch/eval_func/model/proc000000_vp.bin"],
        ["models", "10", "scratch/eval_grad/model/proc000000_vp.bin"],
        ["vtk", "10", "scratch/solver/001/OUTPUT_FILES/DATABASES_MPI/"
                      "proc000000_vp.vtk"]]))


def test_a_file_that_cannot_be_removed_is_reported(workdir):
    found = clean.find_files(os.path.join(workdir, "scratch"),
                             categories=["models"])
    gone = found["models"][0]
    os.remove(gone)
    nfiles, _, failed = clean.remove_files(found, workdir)
    assert((nfiles, failed) == (1, [gone]))


def test_sizes_are_for_a_person():
    assert(clean.human_size(0) == "0 B")
    assert(clean.human_size(999) == "999 B")
    assert(clean.human_size(1500) == "1.5 kB")
    assert(clean.human_size(2.5e9) == "2.5 GB")
    assert(clean.human_size(3e12) == "3000.0 GB")


def test_format_lists_every_file_when_asked(workdir):
    found = clean.find_files(os.path.join(workdir, "scratch"))
    brief = clean.format_files(found, workdir)
    full = clean.format_files(found, workdir, list_files=True)

    assert(f"total: {sum(len(v) for v in LARGE.values())} files" in brief)
    assert("proc000000_save_forward_arrays.bin" not in brief)
    for fid in sum(LARGE.values(), []):
        assert(fid in full)
        assert(fid not in brief)


# -----------------------------------------------------------------------------
# The command
# -----------------------------------------------------------------------------

def _clean_command(workdir, *args, answer=None, capsys=None):
    """Run `seisflows clean` in a working directory, return what it printed"""
    # The paths that `seisflows configure` writes to a parameter file
    with open(os.path.join(workdir, "parameters.yaml"), "w") as f:
        f.write(msg.base_parameter_file)
        for name, path in [("scratch", "scratch"), ("output", "output"),
                           ("log_files", "logs"), ("state_file", "sfstate.txt"),
                           ("output_log", "sflog.txt")]:
            f.write(f"\npath_{name}: {os.path.join(workdir, path)}")
        f.write("\n")
    argv = ["seisflows", "-w", workdir, "-p",
            os.path.join(workdir, "parameters.yaml"), "clean", *args]
    with patch.object(sys, "argv", argv):
        if answer is None:
            SeisFlows()()
        else:
            with patch("builtins.input", return_value=answer) as mock_input:
                SeisFlows()()
            assert(mock_input.called)
    return capsys.readouterr().out


def _exists(workdir, fids):
    return [os.path.exists(os.path.join(workdir, fid)) for fid in fids]


def test_dry_run_removes_nothing_and_names_every_file(workdir, capsys):
    out = _clean_command(workdir, "--dry_run", capsys=capsys)
    for fid in sum(LARGE.values(), []):
        assert(fid in out), fid
        assert(os.path.exists(os.path.join(workdir, fid)))
    for fid in KEPT:
        assert(fid not in out), fid
    assert("nothing was removed" in out)
    assert(not os.path.exists(os.path.join(workdir,
                                           "scratch/system/clean_manifest.txt")))


def test_clean_removes_the_large_files_when_told_to(workdir, capsys):
    out = _clean_command(workdir, answer="y", capsys=capsys)

    assert(not any(_exists(workdir, sum(LARGE.values(), []))))
    assert(all(_exists(workdir, KEPT)))
    assert(all(_exists(workdir, sum(NOT_USUAL.values(), []))))
    assert("removed 20 files" in out)
    assert(os.path.exists(os.path.join(workdir,
                                       "scratch/system/clean_manifest.txt")))


def test_clean_asks_first(workdir, capsys):
    _clean_command(workdir, answer="n", capsys=capsys)
    assert(all(_exists(workdir, sum(LARGE.values(), []))))
    _clean_command(workdir, answer="", capsys=capsys)
    assert(all(_exists(workdir, sum(LARGE.values(), []))))


def test_clean_force_does_not_ask(workdir, capsys):
    _clean_command(workdir, "--force", capsys=capsys)  # no input given
    assert(not any(_exists(workdir, sum(LARGE.values(), []))))


def test_clean_include_and_exclude(workdir, capsys):
    _clean_command(workdir, "-f", "-i", "obs", "-e", "forward_arrays",
                   "kernels", capsys=capsys)

    assert(all(_exists(workdir, LARGE["forward_arrays"] + LARGE["kernels"])))
    assert(not any(_exists(workdir, NOT_USUAL["obs"])))
    assert(not any(_exists(workdir, LARGE["databases"] + LARGE["traces"])))
    assert(all(_exists(workdir, NOT_USUAL["solver_logs"])))


def test_clean_all_is_what_clean_used_to_be(workdir, capsys):
    for name in ["output", "logs"]:
        assert(os.path.isdir(os.path.join(workdir, name)))
    _clean_command(workdir, "--all", answer="y", capsys=capsys)

    assert(os.listdir(workdir) == ["parameters.yaml"])


def test_clean_all_dry_run_removes_nothing(workdir, capsys):
    before = sorted(os.listdir(workdir))
    out = _clean_command(workdir, "--all", "--dry_run", capsys=capsys)

    assert(sorted(os.listdir(workdir)) == before)
    for name in ["scratch", "output", "logs", "sfstate.txt", "sflog.txt"]:
        assert(name in out)
    assert("parameters.yaml" not in out)


def test_clean_says_when_there_is_nothing_to_remove(workdir, capsys):
    _clean_command(workdir, "-f", capsys=capsys)
    out = _clean_command(workdir, "-f", capsys=capsys)
    assert("no files to remove" in out)


def test_clean_needs_a_parameter_file(tmpdir):
    argv = ["seisflows", "-w", str(tmpdir), "-p",
            os.path.join(str(tmpdir), "missing.yaml"), "clean"]
    with patch.object(sys, "argv", argv), pytest.raises(SystemExit):
        SeisFlows()()
