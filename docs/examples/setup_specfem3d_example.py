"""
Script to setup a SPECFEM3D example problem on any system. The only
manual step required by the User is to compile SPECFEM3D binaries.
Should be run directly from the working directory. This script will:

1) Create a new SPECFEM3D working directory 
2) Create a starting model using the Tape 2007 example
3) Create a target model using a modified version of the Tape 2007 example
"""
import os
import glob
import shutil
import subprocess
import numpy as np

from seisflows.tools import unix
from seisflows.seisflows import SeisFlows

# USER CAN EDIT THE FOLLOWING PATHS:
# WORKDIR: points to your own working directory
# SPECFEM3D: points to an existing SPECFEM3D repository if available
WORKDIR = os.getcwd()
SPECFEM3D_ORIGINAL = input("Path to SPECFEM3D respository with compiled "
                           "binaries?: ")

assert(os.path.exists(SPECFEM3D_ORIGINAL)), \
    f"SPECFEM3D repo doesn't exist: {SPECFEM3D_ORIGINAL}"

# Instantiate the SeisFlows class to use its command line arguments
sf = SeisFlows()

# ==============================================================================
# Distribute the necessary file structure of the SPECFEM3D repository that we 
# will downloaded/reference
SPECFEM3D_BIN_ORIGINAL = os.path.join(SPECFEM3D_ORIGINAL, "bin")
SPECFEM3D_DATA_ORIGINAL = os.path.join(SPECFEM3D_ORIGINAL, "DATA")
TAPE_2007_EXAMPLE = os.path.join(SPECFEM3D_ORIGINAL, "EXAMPLES", 
                                 "homogeneous_halfspace")

# The SPECFEM3D working directory that we will create separate from the
# downloaded repo
SPECFEM3D_WORKDIR = os.path.join(WORKDIR, "SPECFEM3D_workdir")
SPECFEM3D_BIN = os.path.join(SPECFEM3D_WORKDIR, "bin")
SPECFEM3D_DATA = os.path.join(SPECFEM3D_WORKDIR, "DATA")
SPECFEM3D_OUTPUT = os.path.join(SPECFEM3D_WORKDIR, "OUTPUT_FILES")

# Pre-defined locations of velocity models we will generate using the solver
SPECFEM3D_MODEL_INIT = os.path.join(SPECFEM3D_WORKDIR, "OUTPUT_FILES_INIT")
SPECFEM3D_MODEL_TRUE = os.path.join(SPECFEM3D_WORKDIR, "OUTPUT_FILES_TRUE")

# ==============================================================================
print("Creating a separate SPECFEM3D working directory")
# Incase we've run this docs page before, delete the working directory before i
# remaking
if os.path.exists(SPECFEM3D_WORKDIR):
    shutil.rmtree(SPECFEM3D_WORKDIR)

os.mkdir(SPECFEM3D_WORKDIR)
os.chdir(SPECFEM3D_WORKDIR)

# Copy the binary files incase we update the source code. These can also be 
# symlinked.
shutil.copytree(SPECFEM3D_BIN_ORIGINAL, "bin")

# Copy the DATA/ directory because we will be making edits here frequently and 
# it's useful to retain the original files for reference. We will be running 
# one of the example problems: Tape2007
shutil.copytree(os.path.join(TAPE_2007_EXAMPLE, "DATA"), "DATA")

# ==============================================================================
if os.path.exists(SPECFEM3D_OUTPUT):
    shutil.rmtree(SPECFEM3D_OUTPUT)

os.mkdir(SPECFEM3D_OUTPUT)

# ==============================================================================
print("Generating initial model from Tape 2007 example")
os.chdir(SPECFEM3D_WORKDIR)

cmd = f"./bin/xdecompose_mesh 4 ./MESH-default {SPECFEM3D_OUTPUT}/DATABASES_MPI"
with open("xdecompose_mesh.out", "w") as f:
    subprocess.run(cmd.split(), shell=True, stdout=f)

cmd = f"mpirun -np 4 ./bin/xgenerate_databases"
with open("xgenerate_databases.out", "w") as f:
    subprocess.run("mpiexec ./bin/xSPECFEM3D", shell=True, stdout=f)

# Move the model files (*.bin) into the OUTPUT_FILES directory, 
# where SeisFlows3 expects them
unix.mv(glob.glob("DATA/*bin"), "OUTPUT_FILES")

# Make sure we don't overwrite this initial model when creating our target model in the next step
unix.mv("OUTPUT_FILES", "OUTPUT_FILES_INIT")

# ==============================================================================
print("Modifying Tape 2007 example for target model")
# GENERATE MODEL_TRUE
os.chdir(SPECFEM3D_DATA)

# Edit the Par_file by increasing velocities by ~10%
sf.sempar("velocity_model", 
          "1 1 2600.d0 5900.d0 3550.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0", 
          skip_print=True)

# ==============================================================================
# Re-run the mesher and solver to generate our target velocity model
os.chdir(SPECFEM3D_WORKDIR)

# Make sure the ./OUTPUT_FILES directory exists since we moved the old one
if os.path.exists(SPECFEM3D_OUTPUT):
    shutil.rmtree(SPECFEM3D_OUTPUT)
os.mkdir(SPECFEM3D_OUTPUT)

print("Generating target model")
# Run the binaries to generate MODEL_TRUE
with open("xmeshfem2D.out", "w") as f:
    subprocess.run("./bin/xmeshfem2D", shell=True, stdout=f)

with open("xSPECFEM3D.out", "w") as f:
    subprocess.run("./bin/xSPECFEM3D", shell=True, stdout=f)

# Move the model files (*.bin) into the OUTPUT_FILES directory, 
# where SeisFlows3 expects them
unix.mv(glob.glob("DATA/*bin"), "OUTPUT_FILES")

# Make sure we don't overwrite this initial model when creating our target model in the next step
unix.mv("OUTPUT_FILES", "OUTPUT_FILES_TRUE")

# ==============================================================================
# Set model==gll in the Par_file so that new models will use .bin files
os.chdir(SPECFEM3D_WORKDIR)
sf.sempar("model", "gll", skip_print=True)

# ==============================================================================
print("Creating SeisFlows3 parameter file")
sf.setup()
