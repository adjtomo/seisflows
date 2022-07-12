"""
Test the ability of the Solver module to interact with various versions of
SPECFEM
"""
import os
import pytest
from glob import glob
from seisflows.tools.specfem import Model
from seisflows.config import ROOT_DIR, NAMES, CFGPATHS


