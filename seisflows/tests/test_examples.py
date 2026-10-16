"""
Test the SPECFEM3D example (Example 4), without needing SPECFEM3D: the
repository is a stand-in with the few files the example reads, and none of the
SPECFEM3D programs are run. See `seisflows examples run 4` for the real thing.
"""
import os
import sys
import pytest
from unittest.mock import patch

from seisflows.tools.config import load_yaml
from seisflows.examples.ex4_hh_specfem3d import SFExample3D

MESH_PAR_FILE = """\
LATITUDE_MIN                    = 0.0
LATITUDE_MAX                    = 134000.0
LONGITUDE_MIN                   = 0.0
LONGITUDE_MAX                   = 134000.0
NEX_XI                          = 36
NEX_ETA                         = 36
NPROC_XI                        = 2
NPROC_ETA                       = 2
# number of materials
NMATERIALS                      = 1
# #material_id  #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id
1   2300.0   2800.0   1500.0   2444.4    300.0 0 2

NREGIONS                        = 1
#NEX_XI_BEGIN  #NEX_XI_END  #NEX_ETA_BEGIN  #NEX_ETA_END  #NZ_BEGIN
1             36             1              36        1        16        1
"""

PAR_FILE = """\
NPROC                           = 4
NSTEP                           = 5000
MODEL                           = default
"""

CMTSOLUTION = """\
PDE  1999 01 01 00 00 00.00  67000 67000 -25000 4.2 4.2 homog_test
event name:       homog_test
time shift:       0.0000
half duration:    5.0
latorUTM:       67000.0
longorUTM:      67000.0
depth:          30.0
Mrr:      -7.600000e+27
"""


@pytest.fixture
def repo(tmpdir):
    """A SPECFEM3D repository with the example and (empty) binaries"""
    repo = os.path.join(str(tmpdir), "specfem3d")
    example = os.path.join(repo, "EXAMPLES", "applications",
                           "homogeneous_halfspace")
    os.makedirs(os.path.join(example, "DATA"))
    os.makedirs(os.path.join(example, "meshfem3D_files"))
    os.makedirs(os.path.join(repo, "bin"))
    for fid, text in [("DATA/Par_file", PAR_FILE),
                      ("DATA/CMTSOLUTION", CMTSOLUTION),
                      ("meshfem3D_files/Mesh_Par_file", MESH_PAR_FILE)]:
        with open(os.path.join(example, fid), "w") as f:
            f.write(text)
    for exc in SFExample3D._binaries:
        open(os.path.join(repo, "bin", exc), "w").close()
    return repo


@pytest.fixture
def workdir(tmpdir, monkeypatch):
    """Run in a clean working directory, the example changes `sys.argv`"""
    monkeypatch.chdir(tmpdir)
    monkeypatch.setattr(sys, "argv", ["seisflows"])
    # The example is meant for a node of a cluster, not for this machine
    monkeypatch.setattr("seisflows.examples.ex4_hh_specfem3d.system_nproc",
                        lambda: 64)
    return str(tmpdir)


def _read(example, *path):
    with open(os.path.join(example.workdir_paths.data, *path), "r") as f:
        return f.read()


@pytest.mark.parametrize("nproc,slices", [
    (1, (1, 1)), (4, (2, 2)), (24, (6, 4)), (28, (7, 4)), (36, (6, 6)),
    (40, (8, 5)), (29, (29, 1)),
])
def test_divide_mesh(nproc, slices):
    """Slices are as square as the number of processors allows"""
    assert(SFExample3D.divide_mesh(nproc) == slices)


@pytest.mark.parametrize("nproc,nex", [(4, (36, 36)), (28, (35, 36)),
                                       (36, (36, 36)), (40, (40, 35))])
def test_mesh_keeps_elements_the_size_of_the_examples(repo, workdir, nproc,
                                                      nex):
    """Every slice has whole elements, which are about 3.7km wide"""
    example = SFExample3D(nproc=nproc, specfem3d_repo=repo)
    assert((example.nex_xi, example.nex_eta) == nex)
    assert(example.nex_xi % example.nproc_xi == 0)
    assert(example.nex_eta % example.nproc_eta == 0)


def test_defaults_fit_a_chinook_node(repo, workdir):
    """28 cores, one event, 25 stations, one iteration, a 5% increase"""
    example = SFExample3D(specfem3d_repo=repo)
    assert((example.nproc, example.ntask, example.nsta, example.niter,
            example.perturb) == (28, 1, 25, 1, 5.))
    assert(example._parameters["components"] == "XYZ")
    assert(example._parameters["nproc"] == 28)


def test_more_processors_than_the_machine_has_is_refused(repo, workdir,
                                                         monkeypatch):
    monkeypatch.setattr("seisflows.examples.ex4_hh_specfem3d.system_nproc",
                        lambda: 8)
    with pytest.raises(AssertionError, match="Run from a compute node"):
        SFExample3D(nproc=28, specfem3d_repo=repo)


def test_repo_can_be_given_with_the_flag_for_specfem2d(repo, workdir):
    """`-r` is the only short flag, so it works for this example too"""
    assert(SFExample3D(specfem2d_repo=repo).sem3d_paths.repo == repo)
    assert(SFExample3D(specfem3d_repo=repo, specfem2d_repo="other"
                       ).sem3d_paths.repo == repo)


def test_missing_binaries_are_listed(repo, workdir, capsys):
    """Nothing is compiled, so say what the User needs to"""
    os.remove(os.path.join(repo, "bin", "xspecfem3D"))
    example = SFExample3D(specfem3d_repo=repo)
    with pytest.raises(SystemExit):
        example.check_specfem3d()
    assert("xspecfem3D" in capsys.readouterr().out)


def test_missing_example_is_reported(tmpdir, workdir, capsys):
    example = SFExample3D(specfem3d_repo=str(tmpdir))
    with pytest.raises(SystemExit):
        example.check_specfem3d()
    assert("--specfem3d_repo" in capsys.readouterr().out)


def test_earlier_workflow_in_the_directory_is_not_continued(repo, workdir,
                                                            capsys):
    """A new parameter file would pick up where the last workflow stopped"""
    open("sfstate.txt", "w").close()
    # Any executable will do for the check that there is MPI
    example = SFExample3D(nproc=4, specfem3d_repo=repo, mpiexec="sh")
    with pytest.raises(SystemExit):
        example.main()
    assert("sfstate.txt" in capsys.readouterr().out)
    assert(not os.path.exists("specfem3d_workdir"))


def test_working_directory_follows_the_arguments(repo, workdir):
    example = SFExample3D(nproc=28, ntask=1, nsta=4, specfem3d_repo=repo)
    example.create_specfem3d_working_directory()

    assert(sorted(os.listdir(example.workdir_paths.bin)) ==
           sorted(SFExample3D._binaries))
    par_file = _read(example, "Par_file")
    assert("NPROC                           = 28" in par_file)
    assert("NSTEP                           = 2000" in par_file)
    mesh = _read(example, "meshfem3D_files", "Mesh_Par_file")
    for line in ["NEX_XI                          = 35",
                 "NEX_ETA                         = 36",
                 "NPROC_XI                        = 7",
                 "NPROC_ETA                       = 4"]:
        assert(line in mesh)
    # Comments are not mistaken for the table of regions
    assert("1  35  1  36  1  16  1\n" in mesh)
    assert("#NEX_XI_BEGIN" in mesh)
    # The material is kept to scale the target model from
    assert(example._material == "1 2300.0 2800.0 1500.0 2444.4 300.0 0 2"
           .split())


def test_target_velocities_are_scaled_and_density_is_not(repo, workdir):
    example = SFExample3D(nproc=4, perturb=-5, specfem3d_repo=repo)
    example.create_specfem3d_working_directory()
    mesh = os.path.join(example.workdir_paths.data, "meshfem3D_files",
                        "Mesh_Par_file")

    values = list(example._material)
    values[2], values[3] = (f"{float(values[2]) * .95:.2f}",
                            f"{float(values[3]) * .95:.2f}")
    example._set_mesh_table_line(mesh, "NMATERIALS", values)
    assert("1  2300.0  2660.00  1425.00  2444.4  300.0  0  2" in
           _read(example, "meshfem3D_files", "Mesh_Par_file"))
    assert(example.perturb == -5.)


def test_events_are_at_the_center_then_on_a_ring(repo, workdir):
    example = SFExample3D(nproc=4, ntask=5, nsta=1, specfem3d_repo=repo)
    example.create_specfem3d_working_directory()

    positions = []
    for i in range(1, 6):
        text = _read(example, f"CMTSOLUTION_{i:0>3}")
        assert(f"homog_test_{i:0>3}" in text)
        assert("depth:          30.0" in text)
        positions.append(tuple(
            float(text.split(key)[1].split()[0])
            for key in ["latorUTM:", "longorUTM:"]))
    assert(positions[0] == (67000., 67000.))
    for lat, lon in positions[1:]:
        assert(((lat - 67000) ** 2 + (lon - 67000) ** 2) ** .5 ==
               pytest.approx(26800., abs=1))
    assert(len(set(positions)) == 5)


@pytest.mark.parametrize("nsta,shape", [(25, (5, 5)), (9, (3, 3)),
                                        (7, (3, 3)), (1, (1, 1))])
def test_stations_are_an_equidistant_grid_inside_the_domain(repo, workdir,
                                                            nsta, shape):
    example = SFExample3D(nproc=4, nsta=nsta, specfem3d_repo=repo)
    example.create_specfem3d_working_directory()

    stations = [line.split() for line in
                _read(example, "STATIONS").strip().splitlines()]
    assert(len(stations) == nsta)
    assert(len({s[0] for s in stations}) == nsta)  # names are unique
    lats = sorted({float(s[2]) for s in stations})
    lons = sorted({float(s[3]) for s in stations})
    assert((len(lats), len(lons)) == shape)
    # Equal spacing, away from the absorbing boundaries on the edges
    for values in [lats, lons]:
        assert(min(values) >= 0.15 * 134000 - 1)
        assert(max(values) <= 0.85 * 134000 + 1)
        steps = [b - a for a, b in zip(values, values[1:])]
        assert(all(step == pytest.approx(steps[0], abs=.1) for step in steps))
    assert(all(s[4:] == ["0.0", "0.0"] for s in stations))  # at the surface


def test_unknown_parameters_are_not_silently_ignored(repo, workdir, capsys):
    """Setting a parameter that doesn't exist does nothing, which hides typos"""
    example = SFExample3D(nproc=4, specfem3d_repo=repo, par=["partiton=x"])
    with pytest.raises(SystemExit):
        example.setup_seisflows_working_directory()
    assert("'partiton' is not a parameter" in capsys.readouterr().out)


def test_user_parameters_override_the_examples(repo, workdir, tmpdir):
    """A YAML file sets many at once, and `par` wins over it"""
    yaml = os.path.join(str(tmpdir), "cluster.yaml")
    with open(yaml, "w") as f:
        f.write("tasktime: 45\ncomponents: Z\nstep_count_max: 3\n")
    example = SFExample3D(nproc=4, specfem3d_repo=repo, par_yaml=yaml,
                          par=["step_count_max=7"])
    example.create_specfem3d_working_directory()
    example.setup_seisflows_working_directory()

    pars = load_yaml("parameters.yaml")
    assert(pars["tasktime"] == 45)
    assert(pars["components"] == "Z")
    assert(pars["step_count_max"] == 7)
    assert(pars["nproc"] == 4)  # the example's, not overridden
    assert(pars["solver"] == "specfem3d")


def test_example_is_listed_and_dispatched(capsys, workdir):
    """`seisflows examples run 4` works because of its file name"""
    from seisflows.seisflows import SeisFlows
    SeisFlows._print_examples()
    assert("4: ex4_hh_specfem3d" in capsys.readouterr().out)


def test_command_line_accepts_the_example_arguments(workdir, monkeypatch):
    from seisflows.seisflows import sfparser
    monkeypatch.setattr(sys, "argv", [
        "seisflows", "examples", "run", "4", "--specfem3d_repo", "/x",
        "--nproc", "28", "--nsta", "25", "--perturb", "-5", "--system",
        "chinook", "--par", "partition=t1small", "walltime=30",
        "--par_yaml", "c.yaml"])
    args = sfparser()[0].parse_args()
    assert((args.choice, args.specfem3d_repo, args.perturb, args.system,
            args.par, args.par_yaml) ==
           (4, "/x", -5., "chinook", ["partition=t1small", "walltime=30"],
            "c.yaml"))
