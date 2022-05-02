#!/usr/bin/env python3
"""
                SEISFLOWS3 SPECFEM2D WORKSTATION EXAMPLE 3

This example will run an inversion between two slightly different homogeneous
halfspace models using 3 sources and 132 stations. It uses the Pyatoa
preprocessing module (which calls the Pyatoa package) to assess misfit and
generate misfit figures.

.. note::
    You can change the number of events (NTASK) and iterations (NITER) by
    changing the constants below the import statements

.. note::
    The tasks involved include:
    1. Download, configure and compile SPECFEM2D
    2. Set up a SPECFEM2D working directory
    3. Generate starting model from Tape2007 example
    4. Generate target model w/ perturbed starting model
    5. Set up a SeisFlows3 working directory
    6. Run two iterations of an inversion workflow

.. rubric::
    $ seisflows examples run 2
"""
import os
import sys
import glob
import shutil
import subprocess

# Re-Use functions from Example 1 to avoid redundant code
from seisflows3.scripts.examples.ex1_specfem2d_workstation_inversion import (
    define_dir_structures, download_specfem2d,
    configure_specfem2d_and_make_binaries, create_specfem2d_working_directory,
    setup_specfem2d_for_model_init, run_xspecfem2d_binaries,
    cleanup_xspecfem2d_run, cleanup_xspecfem2d_run,
    setup_seisflows_working_directory
)
from seisflows3.tools import msg
from seisflows3.config import Dict
from seisflows3.seisflows import SeisFlows
from seisflows3.tools.unix import cd, cp, rm, ln, mv, mkdir

# You are free to change the number of events (NTASK) and iterations (NITER).
# 1 <= NTASK <= 25
# 1 <= NITER <= inf
NTASK = 3
NITER = 2


def setup_specfem2d_for_model_true(sf):
    """
    Make some adjustments to the  parameter file to create the final velocity
    model. This function assumes it is running from inside the SPECFEM2D/DATA
    directory. Also assumes that setup_specfem2d_for_model_init() has been
    run, which set up how the model is output

    :type sf: seisflows3.seisflows.SeisFlows
    :param sf: the SeisFlows3 command line tool
    """
    assert(os.path.exists("Par_file")), f"I cannot find the Par_file!"

    print("> EX3: Updating SPECFEM2D to set checkerboard model as MODEL_TRUE")

    sf.sempar("model", "legacy")  # To read the model_velocity.dat_checker
    rm("proc000000_model_velocity.dat_input")
    ln("model_velocity.dat_checker", "proc000000_model_velocity.dat_input")


def setup_specfem2d_for_seisflows3_inversion(sf):
    """
    Final changes to the SPECFEM2D Par_file before running SeisFlows. Par_file
    will be used to control all the child specfem2d directories. Need to tell
    them to read models from .bin files, and to use existing station files
    rather than create them from the Par_file
    rather than create them from the Par_file
    """
    print("> EX3: Finalizing SPECFEM2D Par_file for SeisFlows3 inversion")
    sf.sempar("model", "gll")  # GLL so SPECFEM reads .bin files
    sf.sempar("use_existing_stations", ".true.")  # Use STATIONS file
    # Assign STATIONS_checker file which has 132 stations
    rm("STATIONS")
    ln("STATIONS_checker", "STATIONS")


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
    sf.par("preprocess", "pyatoa")
    sf.configure()

    sf.par("end", 1)  # only 1 iteration
    sf.par("ntask", ntask)  # we will be using 3 sources for this example
    sf.par("materials", "elastic")  # how the velocity model is parameterized
    sf.par("density", "constant")  # update density or keep constant
    sf.par("nt", 5000)  # set by SPECFEM2D Par_file
    sf.par("dt", .06)  # set by SPECFEM2D Par_file
    sf.par("f0", 0.084)  # set by SOURCE file
    sf.par("format", "ascii")  # how to output synthetic seismograms
    sf.par("case", "synthetic")  # synthetic-synthetic inversion
    sf.par("attenuation", False)
    sf.par("components", "Y")

    # Pyatoa parameters
    sf.par("unit_output", "DISP")

    # The pre-filter here messes with how Pyflex selects windows
    # sf.par("min_period", 1/8)  # slight pre-filter but otherwise data is same
    # sf.par("max_period", 500)


    sf.par("start_pad", 48)  # T0 set in Par_file
    sf.par("end_pad", 5000 * .06)  # nt * dt defined by Par_file

    sf.par("specfem_bin", workdir_paths.bin)
    sf.par("specfem_data", workdir_paths.data)
    sf.par("model_init", workdir_paths.model_init)
    sf.par("model_true", workdir_paths.model_true)


def main(run_example=False):
    """
    Setup the example and then optionally run the actual seisflows workflow

    TODO Figure out how to not have to repeat this code block, which is
    TODO an exact copy of example 1

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
    setup_specfem2d_for_seisflows3_inversion(sf3_cli_tool)
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
                  f"between a homogeneous halfspace model and a checkerboard "
                  f"model using {NTASK} sources and 132 stations."
                  "The tasks involved include: ",
                  items=["1. (optional) Download, configure, compile SPECFEM2D",
                         "2. Set up a SPECFEM2D working directory",
                         "3. Generate starting model from Tape2007 example",
                         "4. Generate target checkerboard model",
                         "5. Set up a SeisFlows3 working directory",
                         f"6. Run {NITER} iterations of an inversion workflow"],
                  header="seisflows3 example 2",
                  border="="))

    if len(sys.argv) > 1:
        main(run_example=sys.argv[1])
