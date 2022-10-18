"""
Test the system modules ability to run jobs
"""
import os
import pytest
from glob import glob
from seisflows.tools import unix
from seisflows.tools.config import get_task_id
from seisflows.system.workstation import Workstation
from seisflows.system.cluster import Cluster


def _test_function_a(tmpdir):
    """
    A test function that simply prints, used to check run
    """
    test_dir = os.path.join(tmpdir, "test_dir")
    if not os.path.exists(test_dir):
        unix.mkdir(test_dir)


def _test_function_b(tmpdir):
    """
    A test function that simply prints
    """
    idx = get_task_id()
    test_fid = f"test_file_{idx:0>2}.txt"
    test_file = os.path.join(tmpdir, "test_dir", test_fid)
    print(f"making test file {test_file}")
    open(test_file, "w").close()


def test_workstation_run(tmpdir, ntask=26):
    """
    Make sure that Workstation run simply calls the functions in question
    """
    system = Workstation(workdir=tmpdir, path_system=tmpdir, ntask=ntask)
    system.setup()  # make required dir. structure
    system.run(funcs=[_test_function_a, _test_function_b], tmpdir=tmpdir)

    assert(len(glob(os.path.join(system.path.log_files, "*"))) == ntask)
    assert(len(glob(os.path.join(tmpdir, "test_dir", "*"))) == ntask)


def test_workstation_run_single(tmpdir, ntask=26):
    """
    Make sure that Workstation run simply calls the functions only once
    """
    system = Workstation(workdir=tmpdir, path_system=tmpdir, ntask=ntask)
    system.setup()  # make required dir. structure
    system.run(funcs=[_test_function_a, _test_function_b], tmpdir=tmpdir,
               single=True)

    assert(len(glob(os.path.join(system.path.log_files, "*"))) == 1)
    assert(len(glob(os.path.join(tmpdir, "test_dir", "*"))) == 1)


def test_cluster_run(tmpdir, ntask=26):
    """
    Make sure that Cluster run can pickle the tasks and run them with subprocess
    """
    system = Cluster(workdir=tmpdir, path_system=tmpdir, ntask=ntask)
    system.setup()  # make required dir. structure
    system.run(funcs=[_test_function_a, _test_function_b], tmpdir=tmpdir)

    assert(len(glob(os.path.join(system.path.log_files, "*"))) == ntask)
    assert(len(glob(os.path.join(tmpdir, "test_dir", "*"))) == ntask)


def test_cluster_run_single(tmpdir, ntask=26):
    """
    Make sure that Cluster runs single tasks accetpably
    """
    system = Cluster(workdir=tmpdir, path_system=tmpdir, ntask=ntask)
    system.setup()  # make required dir. structure
    system.run(funcs=[_test_function_a, _test_function_b], tmpdir=tmpdir,
               single=True)

    assert(len(glob(os.path.join(system.path.log_files, "*"))) == 1)
    assert(len(glob(os.path.join(tmpdir, "test_dir", "*"))) == 1)


