"""
Create a test directory structure with actual data to test individual modules
"""
import os
from seisflows.tools import unix
from seisflows.tools.config import CFGPATHS, ROOT_DIR


testdir = os.path.join(ROOT_DIR, "tests", "test_data")
workdir = os.path.join(testdir, "workdir")
for key, path in CFGPATHS.__dict__.items():
    full_path = os.path.join(workdir, path)
    if os.path.exists(full_path):
        unix.rm(full_path)
    # Identify files by the separator between filename and extension
    if "." in path:
        open(full_path, "w")
    else:
        unix.mkdir(os.path.join(workdir, path))

# Make the scratch solver directory
scratchdir = os.path.join(workdir, CFGPATHS.SCRATCHDIR)
solverdir = os.path.join(scratchdir, "solver", "001")
tracesdir = os.path.join(solverdir, "traces")
for dir_ in [scratchdir, solverdir, tracesdir]:
    unix.mkdir(dir_)

for dir_ in ["obs", "syn", "adj"]:
    unix.mkdir(os.path.join(tracesdir, dir_))
