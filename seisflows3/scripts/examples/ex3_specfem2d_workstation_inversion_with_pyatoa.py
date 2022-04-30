#!/usr/bin/env python3
"""
                SEISFLOWS3 SPECFEM2D WORKSTATION EXAMPLE 1

This example will run two iterations of an inversion to assess misfit between
a homogeneous halfspace model and a slightly perturbed homogeneous halfspace
model using 3 events and 1 receiver.

.. note::
    You can change the number of events (NTASK) and iterations (NITER) by
    changing the constants below the import statements

.. warning::
    Because we are using 3 events and only 1 receiver, the results of this
    inversion will be of questionable quality. This example is only meant
    to highlight how SeisFlows3 operates during an inversion workflow.

.. note::
    The tasks involved include:
    1. Download, configure and compile SPECFEM2D
    2. Set up a SPECFEM2D working directory
    3. Generate starting model from Tape2007 example
    4. Generate target model w/ perturbed starting model
    5. Set up a SeisFlows3 working directory
    6. Run two iterations of an inversion workflow

.. rubric::
    $ python ex1_specfem2d_workstation_inversion.py

    OR

    $ seisflows examples run 1
"""
import os
import sys
import glob
import shutil
import subprocess

from seisflows3.tools import msg
from seisflows3.config import Dict
from seisflows3.seisflows import SeisFlows
from seisflows3.tools.unix import cd, cp, rm, ln, mv, mkdir

# You are free to change the number of events (NTASK) and iterations (NITER).
# 1 <= NTASK <= 25
# 1 <= NITER <= inf
NTASK = 3
NITER = 2


def define_dir_structures(cwd, specfem2d_repo, ex="Tape2007"):
    """
    Define the example directory structure, which will contain abridged
    versions of the SPECFEM2D working directory

    :type cwd: str
    :param cwd: current working directory
    :type specfem2d_repo: str
    :param specfem2d_repo: location of the SPECFEM2D repository
    :type ex: str
    :type ex: The name of the example problem inside SPECFEM2D/EXAMPLES
    """
    if not specfem2d_repo:
        print(f"No existing SPECFEM2D repo given, default to: {cwd}/specfem2d")
        specfem2d_repo = os.path.join(cwd, "specfem2d")
    else:
        assert(os.path.exists(specfem2d_repo)), (
            f"User supplied SPECFEM2D directory '{specfem2d_repo}' "
            f"does not exist, please check your path and try again."
        )

    # This defines required structures from the SPECFEM2D repository
    sem2d = {
        "repo": specfem2d_repo,
        "bin": os.path.join(specfem2d_repo, "bin"),
        "data": os.path.join(specfem2d_repo, "DATA"),
        "example": os.path.join(specfem2d_repo, "EXAMPLES", ex),
        "example_data": os.path.join(specfem2d_repo, "EXAMPLES", ex, "DATA")
    }
    assert(os.path.exists(sem2d["example"])), (
        f"SPECFEM2D/EXAMPLE directory: '{sem2d['example']}' "
        f"does not exist, please check this path and try again."
    )
    # This defines a working directory structure which we will create
    working_directory = os.path.join(cwd, "specfem2d_workdir")
    workdir = {
        "workdir": os.path.join(working_directory),
        "bin": os.path.join(working_directory, "bin"),
        "data": os.path.join(working_directory, "DATA"),
        "output": os.path.join(working_directory, "OUTPUT_FILES"),
        "model_init": os.path.join(working_directory, "OUTPUT_FILES_INIT"),
        "model_true": os.path.join(working_directory, "OUTPUT_FILES_TRUE"),
    }
    return Dict(sem2d), Dict(workdir)


def download_specfem2d():
    """
    Download the latest version of SPECFEM2D from GitHub, devel branch.
    Last successfully tested 4/28/22
    """
    cmd = ("git clone --recursive --branch devel " 
           "https://github.com/geodynamics/specfem2d.git")

    print(f"Downloading SPECFEM2D with command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def configure_specfem2d_and_make_binaries():
    """
    Run ./configure within the SPECFEM2D repo directory.
    This function assumes it is being run from inside the repo. Should guess
    all the configuration options. Probably the least stable part of the example
    """
    if not os.path.exists("./config.log"):
        cmd = "./configure"
        print(f"Configuring SPECFEM2D with command: {cmd}")
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
    else:
        print("SPECFEM2D already configured, skipping 'configure'")

    if not glob.glob("./bin/x*"):
        cmd = "make all"
        print(f"Making SPECFEM2D binaries with command: {cmd}")
        subprocess.run(cmd, shell=True, check=True,  stdout=subprocess.DEVNULL)
    else:
        print("executables found in SPECFEM2D/bin directory, skipping 'make'")


def create_specfem2d_working_directory(sem2d_paths, workdir_paths):
    """
    Create the working directory where we will generate our initial and
    final models using one of the SPECFEM2D examples

    :type sem2d_paths: Dict
    :param sem2d_paths: path dictionary for the SPECFEM2D repository
    :type workdir_paths: Dict
    :param workdir_paths: path dictionary for the SPECFEM2D working directory
    """
    # Incase this example has written before, remove dir. that were created
    rm(workdir_paths.workdir)
    mkdir(workdir_paths.workdir)

    # Copy the binary executables and DATA from the SPECFEM2D example
    cp(sem2d_paths.bin, workdir_paths.bin)
    cp(sem2d_paths.example_data, workdir_paths.data)

    # Make sure that SPECFEM2D can find the expected files in the DATA/ dir
    cd(workdir_paths.data)
    rm("SOURCE")
    ln("SOURCE_001", "SOURCE")
    rm("Par_file")
    ln("Par_file_Tape2007_onerec", "Par_file")


def setup_specfem2d_for_model_init(sf):
    """
    Make some adjustments to the original parameter file to.
    This function assumes it is running from inside the SPECFEM2D/DATA directory

    :type sf: seisflows3.seisflows.SeisFlows
    :param sf: the SeisFlows3 command line tool
    """
    assert(os.path.exists("Par_file")), f"I cannot find the Par_file!"

    print("> Setting the SPECFEM2D Par_file for SeisFlows3 compatiblility")

    sf.sempar("setup_with_binary_database", 1)  # allow creation of .bin files
    sf.sempar("save_model", "binary")  # output model in .bin database format
    sf.sempar("save_ASCII_kernels", ".false.")  # output kernels in .bin format


def setup_specfem2d_for_model_true(sf):
    """
    Make some adjustments to the  parameter file to create the final velocity
    model. This function assumes it is running from inside the SPECFEM2D/DATA
    directory

    :type sf: seisflows3.seisflows.SeisFlows
    :param sf: the SeisFlows3 command line tool
    """
    assert(os.path.exists("Par_file")), f"I cannot find the Par_file!"

    print("> Updating initial homogeneous velocity model values")

    new_model = "1 1 2600.d0 5900.d0 3550.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0"

    sf.sempar("velocity_model", new_model)


def run_xspecfem2d_binaries():
    """
    Runs the xmeshfem2d and then xspecfem2d binaries using subprocess and then
    do some cleanup to get files in the correct locations. Assumes that we
    can run the binaries directly with './'
    """
    cmd_mesh = f"./bin/xmeshfem2D > OUTPUT_FILES/mesher.log.txt"
    cmd_spec = f"./bin/xspecfem2D > OUTPUT_FILES/solver.log.txt"

    for cmd in [cmd_mesh, cmd_spec]:
        print(f"Running SPECFEM2D with command: {cmd}")
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)


def cleanup_xspecfem2d_run(workdir_paths, new_name=None):
    """
    Do some cleanup after running the SPECFEM2D binaries to make sure files are
    in the correct locations, and rename the OUTPUT_FILES directory so that it
    does not get overwritten by subsequent runs

    :type workdir_paths: Dict
    :param workdir_paths: path dictionary for the SPECFEM2D working directory
    :type new_name: str
    :param new_name: if we want to rename OUTPUT_FILES to something else,
        so that it does not get overwritten by subsequent simulations
    """
    print("> Cleaning up after xspecfem2d, moving files to correct locations")
    cd(workdir_paths.workdir)
    mv(glob.glob("DATA/*bin"), workdir_paths.output)

    if new_name is not None:
        mv(workdir_paths.output, new_name)


def setup_seisflows_working_directory(sf, workdir_paths, ntask=3, niter=2):
    """
    Create and set the SeisFlows3 parameter file, making sure all required
    parameters are set correctly for this example problem

    :type sf: seisflows3.seisflows.SeisFlows
    :param sf: the SeisFlows3 command line tool
    :type workdir_paths: Dict
    :param workdir_paths: path dictionary for the SPECFEM2D working directory
    :type ntask: int
    :param ntask: Number of sources to include in the inversion
    :type niter: int
    :param ninter: number of iterations to perform within the inversion
    """
    sf.setup(force=True)  # Force will delete any existing parameter file
    sf.configure()

    sf.par("ntask", ntask)  # we will be using 3 sources for this example
    sf.par("materials", "elastic")  # how the velocity model is parameterized
    sf.par("density", "constant")  # update density or keep constant
    sf.par("nt", 5000)  # set by SPECFEM2D Par_file
    sf.par("dt", .06)  # set by SPECFEM2D Par_file
    sf.par("f0", 0.084)  # set by SOURCE file
    sf.par("format", "ascii")  # how to output synthetic seismograms
    sf.par("begin", 1)  # first iteration
    sf.par("end", niter)  # final iteration -- we will run 2
    sf.par("case", "synthetic")  # synthetic-synthetic inversion
    sf.par("attenuation", False)

    sf.par("specfem_bin", workdir_paths.bin)
    sf.par("specfem_data", workdir_paths.data)
    sf.par("model_init", workdir_paths.model_init)
    sf.par("model_true", workdir_paths.model_true)


def main(run_example=False):
    """
    Setup the example and then optionally run the actual seisflows workflow

    :type run_example: bool
    :param run_example: Directly run the seisflows workflow after the setup
    """
    sys.argv = [sys.argv[0]]  # Ensure that no arguments are given to the CLI
    sf3_cli_tool = SeisFlows()
    cwd = os.getcwd()
    specfem2d_repo = input(
        msg.cli("If you have already downloaded SPECMFE2D, please input its "
                "path here. If blank, this example will pull the latest version "
                "from GitHub and attempt to configure and make the "
                "binaries:\n> ")
    )

    # Step 0: Initialize path structure and the SeisFlows command line tool
    print(msg.cli("EXAMPLE SETUP", border="="))
    sem2d_paths, workdir_paths = define_dir_structures(cwd, specfem2d_repo)

    # Step 1: Download and configure SPECFEM2D, make binaries. Or check if done
    if not os.path.exists(sem2d_paths.repo):
        download_specfem2d()

    cd(sem2d_paths.repo)
    try:
        configure_specfem2d_and_make_binaries()
    except subprocess.CalledProcessError as e:
        print(f"configure and make step has failed, please check and retry. "
              f"If this command repeatedly fails, you may need to "
              f"configure and compile SPECFEM2D manually.\n{e}")
        sys.exit(-1)

    # Step 2: Create a working directory and generate the initial/final models
    print(msg.cli("SETTING UP SPECFEM2D WORKING DIRECTORY", border="="))
    create_specfem2d_working_directory(sem2d_paths, workdir_paths)

    # Step 2a: Par_file manipulations and prepare directory structure
    cd(workdir_paths.data)
    setup_specfem2d_for_model_init(sf3_cli_tool)
    rm(workdir_paths.output)
    mkdir(workdir_paths.output)

    # Step 2b: Generate MODEL_INIT, rearrange consequent directory structure
    cd(workdir_paths.workdir)
    print(msg.cli("GENERATING INITIAL MODEL", border="="))
    run_xspecfem2d_binaries()
    cleanup_xspecfem2d_run(workdir_paths, new_name=workdir_paths.model_init)

    # Step 3b: Prepare Par_file and directory for MODEL_TRUE generation
    cd(workdir_paths.data)
    setup_specfem2d_for_model_true(sf3_cli_tool)
    rm(workdir_paths.output)
    mkdir(workdir_paths.output)

    # Step 3d: Generate MODEL_TRUE, rearrange consequent directory structure
    print(msg.cli("GENERATING TRUE MODEL", border="="))
    cd(workdir_paths.workdir)
    run_xspecfem2d_binaries()
    cleanup_xspecfem2d_run(workdir_paths, new_name=workdir_paths.model_true)

    # Step 3: Setup and run SeisFlows3 to perform 2 iterations of an inversion
    #   Using the velocity models we just generated
    print(msg.cli("SETTING UP SEISFLOWS3", border="="))
    cd(cwd)
    setup_seisflows_working_directory(sf3_cli_tool, workdir_paths, ntask=NTASK,
                                      niter=NITER)

    # Step 3b: last minute Par_file edits
    cd(workdir_paths.data)
    sf3_cli_tool.sempar("model", "gll")  # GLL so SPECFEM reads .bin files
    print("> setup complete")

    if run_example == "run":
        # Step 4: Run the inversion!
        print(msg.cli("RUNNING SEISFLOWS3 INVERSION WORKFLOW", border="="))
        # Run SeisFlows3 Inversion as an external process so that it doesn't get
        # affected by our current environment
        cd(cwd)
        subprocess.run("seisflows submit -f", check=False, shell=True)


if __name__ == "__main__":
    # Some header information before starting to inform the user of goings ons
    print(msg.ascii_logo_small)
    print(msg.cli(f"This is a [SPECFEM2D] [WORKSTATION] example, which will "
                  f"run {NITER} iterations of an inversion to assess misfit "
                  f"between two homogeneous halfspace models with slightly "
                  f"different velocities, {NTASK} sources and 1 receiver. "
                  f"The tasks involved include: ",
                  items=["1. (optional) Download, configure, compile SPECFEM2D",
                         "2. Set up a SPECFEM2D working directory",
                         "3. Generate starting model from Tape2007 example",
                         "4. Generate target model w/ perturbed starting model",
                         "5. Set up a SeisFlows3 working directory",
                         f"6. Run {NITER} iterations of an inversion workflow"],
                  header="seisflows3 example 1",
                  border="="))

    if len(sys.argv) > 1:
        main(run_example=sys.argv[1])

