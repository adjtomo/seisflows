"""
Test the ability of the Solver module to interact with various versions of
SPECFEM
"""
import os
import pytest
from glob import glob
from seisflows.tools.specfem import Model
from seisflows.config import ROOT_DIR, NAMES, CFGPATHS

from seisflows.solver.specfem import Specfem
from seisflows.solver.specfem2d import Specfem2D
from seisflows.solver.specfem3d import Specfem3D
from seisflows.solver.specfem3d_globe import Specfem3DGlobe


TEST_DATA = os.path.join(ROOT_DIR, "tests", "test_data", "test_solver")


def test_taskid():
    """
    Make sure that task id returns correctly
    """
    solver = Specfem()
    assert(solver.taskid == 0)
    os.environ["SEISFLOWS_TASKID"] = "9"
    assert(solver.taskid == 9)


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


def test_data_filenames():
    """
    Test that data filenames are returned correctly
    """
    sources = os.path.join(TEST_DATA, "test_solver", "sources")
    solver = Specfem(path_specfem_data=sources, source_prefix="CMTSOLUTION")
    pytest.set_trace()