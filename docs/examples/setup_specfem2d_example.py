"""
Script to setup a Specfem2D example problem on any system. The only
manual step required by the User is to compile Specfem2D binaries.
Should be run directly from the working directory. This script will:

1) Create a new SPECFEM2D working directory 
2) Create a starting model using the Tape 2007 example
3) Create a target model using a modified version of the Tape 2007 example
"""
import os
import glob
import shutil
import subprocess

from seisflows.tools import unix
from seisflows.seisflows import SeisFlows

# USER CAN EDIT THE FOLLOWING PATHS:
# WORKDIR: points to your own working directory
# SPECFEM2D: points to an existing specfem2D repository if available
WORKDIR = os.getcwd()
SPECFEM2D_ORIGINAL = input("Path to SPECFEM2D respository with compiled "
                           "binaries?: ")

assert(os.path.exists(SPECFEM2D_ORIGINAL)), \
    f"SPECFEM2D repo doesn't exist: {SPECFEM2D_ORIGINAL}"

# Instantiate the SeisFlows class to use its command line arguments
sf = SeisFlows()

# ==============================================================================
# Distribute the necessary file structure of the SPECFEM2D repository that we 
# will downloaded/reference
SPECFEM2D_BIN_ORIGINAL = os.path.join(SPECFEM2D_ORIGINAL, "bin")
SPECFEM2D_DATA_ORIGINAL = os.path.join(SPECFEM2D_ORIGINAL, "DATA")
TAPE_2007_EXAMPLE = os.path.join(SPECFEM2D_ORIGINAL, "EXAMPLES", "Tape2007")

# The SPECFEM2D working directory that we will create separate from the
# downloaded repo
SPECFEM2D_WORKDIR = os.path.join(WORKDIR, "specfem2d_workdir")
SPECFEM2D_BIN = os.path.join(SPECFEM2D_WORKDIR, "bin")
SPECFEM2D_DATA = os.path.join(SPECFEM2D_WORKDIR, "DATA")
SPECFEM2D_OUTPUT = os.path.join(SPECFEM2D_WORKDIR, "OUTPUT_FILES")

# Pre-defined locations of velocity models we will generate using the solver
SPECFEM2D_MODEL_INIT = os.path.join(SPECFEM2D_WORKDIR, "OUTPUT_FILES_INIT")
SPECFEM2D_MODEL_TRUE = os.path.join(SPECFEM2D_WORKDIR, "OUTPUT_FILES_TRUE")

# ==============================================================================
print("Creating a separate SPECFEM2D working directory")
# Incase we've run this docs page before, delete the working directory before i
# remaking
if os.path.exists(SPECFEM2D_WORKDIR):
    shutil.rmtree(SPECFEM2D_WORKDIR)

os.mkdir(SPECFEM2D_WORKDIR)
os.chdir(SPECFEM2D_WORKDIR)

# Copy the binary files incase we update the source code. These can also be 
# symlinked.
shutil.copytree(SPECFEM2D_BIN_ORIGINAL, "bin")

# Copy the DATA/ directory because we will be making edits here frequently and 
# it's useful to retain the original files for reference. We will be running 
# one of the example problems: Tape2007
shutil.copytree(os.path.join(TAPE_2007_EXAMPLE, "DATA"), "DATA")

# ==============================================================================
# First we will set the correct SOURCE and STATION files.
# This is the same task as shown in ./run_this_example.sh
os.chdir(SPECFEM2D_DATA)

# Symlink source 001 as our main source
if os.path.exists("SOURCE"):
    os.remove("SOURCE")
os.symlink("SOURCE_001", "SOURCE")

# Copy the correct Par_file so that edits do not affect the original file
if os.path.exists("Par_file"):
    os.remove("Par_file")
shutil.copy("Par_file_Tape2007_onerec", "Par_file")

# ==============================================================================
print("Adjusting SPECFEM2D parameter file")
# Ensure that parameter files are set correctly
os.chdir(SPECFEM2D_DATA)

sf.sempar("setup_with_binary_database", 1, skip_print=True)
sf.sempar("save_model", "binary", skip_print=True)
sf.sempar("save_ascii_kernels", ".false.", skip_print=True)

os.chdir(SPECFEM2D_WORKDIR)

if os.path.exists(SPECFEM2D_OUTPUT):
    shutil.rmtree(SPECFEM2D_OUTPUT)

os.mkdir(SPECFEM2D_OUTPUT)

# ==============================================================================
print("Generating initial model from Tape 2007 example")
os.chdir(SPECFEM2D_WORKDIR)

with open("xmeshfem2D.out", "w") as f:
    subprocess.run("./bin/xmeshfem2D", shell=True, stdout=f)

with open("xspecfem2D.out", "w") as f:
    subprocess.run("./bin/xspecfem2D", shell=True, stdout=f)

# Move the model files (*.bin) into the OUTPUT_FILES directory, 
# where SeisFlows3 expects them
unix.mv(glob.glob("DATA/*bin"), "OUTPUT_FILES")

# Make sure we don't overwrite this initial model when creating our target
# model in the next step
unix.mv("OUTPUT_FILES", "OUTPUT_FILES_INIT")

# ==============================================================================
print("Modifying Tape 2007 example for target model")
# GENERATE MODEL_TRUE
os.chdir(SPECFEM2D_DATA)

# Edit the Par_file by increasing velocities by ~10%
sf.sempar("velocity_model", 
          "1 1 2600.d0 5900.d0 3550.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0", 
          skip_print=True)

# ==============================================================================
# Re-run the mesher and solver to generate our target velocity model
os.chdir(SPECFEM2D_WORKDIR)

# Make sure the ./OUTPUT_FILES directory exists since we moved the old one
if os.path.exists(SPECFEM2D_OUTPUT):
    shutil.rmtree(SPECFEM2D_OUTPUT)
os.mkdir(SPECFEM2D_OUTPUT)

print("Generating target model")
# Run the binaries to generate MODEL_TRUE
with open("xmeshfem2D.out", "w") as f:
    subprocess.run("./bin/xmeshfem2D", shell=True, stdout=f)

with open("xspecfem2D.out", "w") as f:
    subprocess.run("./bin/xspecfem2D", shell=True, stdout=f)

# Move the model files (*.bin) into the OUTPUT_FILES directory, 
# where SeisFlows3 expects them
unix.mv(glob.glob("DATA/*bin"), "OUTPUT_FILES")
unix.mv("OUTPUT_FILES", "OUTPUT_FILES_TRUE")

# ==============================================================================
# Set model==gll in the Par_file so that new models will use .bin files
os.chdir(SPECFEM2D_WORKDIR)
sf.sempar("model", "gll", skip_print=True)

# ==============================================================================
print("Creating SeisFlows3 parameter file")
sf.setup()
