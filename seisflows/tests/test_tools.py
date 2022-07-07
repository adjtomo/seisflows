"""
Test any of the utility functions defined in the Tools directory
"""
import os
import pytest
from seisflows.config import ROOT_DIR, NAMES, CFGPATHS


TEST_DIR = os.path.join(ROOT_DIR, "tests")


def test_load_specfem_model():
    """
    Make sure we can dynamically load SPECFEM models in various formats
    """
