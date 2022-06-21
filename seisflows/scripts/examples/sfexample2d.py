#!/usr/bin/env python3
"""
                SEISFLOWS SPECFEM2D WORKSTATION EXAMPLE 1

This example will run two iterations of an inversion to assess misfit between
a homogeneous halfspace model and a slightly perturbed homogeneous halfspace
model using 3 events and 1 receiver.

.. note::
    You can change the number of events (NTASK) and iterations (NITER) by
    changing the constants below the import statements

.. warning::
    Because we are using 3 events and only 1 receiver, the results of this
    inversion will be of questionable quality. This example is only meant
    to highlight how SeisFlows operates during an inversion workflow.

.. note::
    The tasks involved include:
    1. Download, configure and compile SPECFEM2D
    2. Set up a SPECFEM2D working directory
    3. Generate starting model from Tape2007 example
    4. Generate target model w/ perturbed starting model
    5. Set up a SeisFlows working directory
    6. Run two iterations of an inversion workflow

.. rubric::
    $ seisflows examples run 1
"""
import os
import sys
import glob
import shutil
import subprocess
import numpy as np

from seisflows.tools import msg
from seisflows.config import Dict
from seisflows.seisflows import SeisFlows
from seisflows.tools.unix import cd, cp, rm, ln, mv, mkdir


class SFExample2D:
    """
    A class for running SeisFlows examples. Simplifies calls structure so that
    multiple example runs can benefit from the code written here
    """
    def __init__(self, ntask=3, niter=2):
        """
        Set path structure which is used to navigate around SPECFEM repositories
        and the example working directory

        :type ntask: int
        :param ntask: number of events to use in inversion, between 1 and 25.
            defaults to 3
        :type niter: int
        :param niter: number of iterations to run. defaults to 2
        """
        specfem2d_repo = input(
            msg.cli("If you have already downloaded SPECMFE2D, please input "
                    "the full path to the repo. If left blank, this example "
                    "will pull the latest version from GitHub and attempt "
                    "to configure and make the binaries:\n> ")
        )

        self.cwd = os.getcwd()
        self.sem2d_paths, self.workdir_paths = self.define_dir_structures(
            cwd=self.cwd, specfem2d_repo=specfem2d_repo
        )
        self.ntask = ntask
        assert(1 <= self.ntask <= 25), \
            f"number of tasks/events must be between 1 and 25, not {self.ntask}"
        self.niter = niter
        assert(1 <= self.niter <= np.inf), \
            f"number of iterations must be between 1 and inf, not {self.niter}"

        # This bool information is provided by the User running 'setup' or 'run'
        self.run_example = bool(sys.argv[1] == "run")

        # Command line tool to use $ seisflows <cmd> from inside Python
        # Zero out sys.argv to ensure that no arguments are given to the CLI
        sys.argv = [sys.argv[0]]
        self.sf = SeisFlows()

    @staticmethod
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
            print(f"No existing SPECFEM2D repo given, default to: "
                  f"{cwd}/specfem2d")
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

    def download_specfem2d(self):
        """
        Download the latest version of SPECFEM2D from GitHub, devel branch.
        Last successfully tested 4/28/22
        """
        if not os.path.exists(self.sem2d_paths.repo):
            cmd = ("git clone --recursive --branch devel " 
                   "https://github.com/geodynamics/specfem2d.git")

            print(f"Downloading SPECFEM2D with command: {cmd}")
            subprocess.run(cmd, shell=True, check=True)

    def configure_specfem2d_and_make_binaries(self):
        """
        Run ./configure within the SPECFEM2D repo directory.
        This function assumes it is being run from inside the repo. Should guess
        all the configuration options. Probably the least stable part of the
        example
        """
        cd(self.sem2d_paths.repo)
        try:
            if not os.path.exists("./config.log"):
                cmd = "./configure"
                print(f"Configuring SPECFEM2D with command: {cmd}")
                # Ignore the configure outputs from SPECFEM
                subprocess.run(cmd, shell=True, check=True,
                               stdout=subprocess.DEVNULL)
            else:
                print("SPECFEM2D already configured, skipping 'configure'")
        except subprocess.CalledProcessError as e:
            print(f"SPECFEM `configure` step has failed, please check and "
                  f"retry. If this command repeatedly fails, you may need "
                  f"to configure SPECFEM2D manually.\n{e}")
            sys.exit(-1)

        try:
            if not glob.glob("./bin/x*"):
                cmd = "make all"
                print(f"Making SPECFEM2D binaries with command: {cmd}")
                # Ignore the make outputs from SPECFEM
                subprocess.run(cmd, shell=True, check=True,
                               stdout=subprocess.DEVNULL)
            else:
                print("executables found in SPECFEM2D/bin directory, "
                      "skipping 'make'")
        except subprocess.CalledProcessError as e:
            print(f"SPECFEM 'make' step has failed, please check and "
                  f"retry. If this command repeatedly fails, you may need "
                  f"to compile SPECFEM2D manually.\n{e}")
            sys.exit(-1)

        # Symlink the finished repo into the working directory so that any
        # subsequent runs won't need to have the user re-type repo location
        if not os.path.exists(os.path.join(self.cwd, "specfem2d")):
            print("symlinking existing specfem2D repository to cwd")
            ln(self.sem2d_paths.repo, os.path.join(self.cwd, "specfem2d"))

    def create_specfem2d_working_directory(self):
        """
        Create the working directory where we will generate our initial and
        final models using one of the SPECFEM2D examples
        """
        assert(os.path.exists(self.sem2d_paths["example"])), (
            f"SPECFEM2D/EXAMPLE directory: '{self.sem2d['example']}' "
            f"does not exist, please check this path and try again."
        )

        # Incase this example has written before, remove dir. that were created
        rm(self.workdir_paths.workdir)
        mkdir(self.workdir_paths.workdir)

        # Copy the binary executables and DATA from the SPECFEM2D example
        cp(self.sem2d_paths.bin, self.workdir_paths.bin)
        cp(self.sem2d_paths.example_data, self.workdir_paths.data)

        # Make sure that SPECFEM2D can find the expected files in the DATA/ dir
        cd(self.workdir_paths.data)
        rm("SOURCE")
        ln("SOURCE_001", "SOURCE")
        rm("Par_file")
        ln("Par_file_Tape2007_onerec", "Par_file")

    def setup_specfem2d_for_model_init(self):
        """
        Make some adjustments to the original parameter file to.
        This function assumes it is running from inside the SPECFEM2D/DATA dir

        """
        cd(self.workdir_paths.data)
        assert(os.path.exists("Par_file")), f"I cannot find the Par_file!"

        print("> Setting the SPECFEM2D Par_file for SeisFlows compatiblility")

        self.sf.sempar("setup_with_binary_database", 1)  # create .bin files
        self.sf.sempar("save_model", "binary")  # output model in .bin format
        self.sf.sempar("save_ASCII_kernels", ".false.")  # kernels also .bin

        rm(self.workdir_paths.output)
        mkdir(self.workdir_paths.output)

    def setup_specfem2d_for_model_true(self):
        """
        Make some adjustments to the parameter file to create the final velocity
        model. This function assumes it is running from inside the
        SPECFEM2D/DATA directory
        """
        cd(self.workdir_paths.data)
        assert(os.path.exists("Par_file")), f"I cannot find the Par_file!"

        print("> Updating initial homogeneous velocity model values")
        new_model = "1 1 2600.d0 5900.d0 3550.0d0 0 0 10.d0 10.d0 0 0 0 0 0 0"
        self.sf.sempar("velocity_model", new_model)

    def run_xspecfem2d_binaries(self):
        """
        Runs the xmeshfem2d and then xspecfem2d binaries using subprocess and then
        do some cleanup to get files in the correct locations. Assumes that we
        can run the binaries directly with './'
        """
        cd(self.workdir_paths.workdir)

        cmd_mesh = f"./bin/xmeshfem2D > OUTPUT_FILES/mesher.log.txt"
        cmd_spec = f"./bin/xspecfem2D > OUTPUT_FILES/solver.log.txt"

        for cmd in [cmd_mesh, cmd_spec]:
            print(f"Running SPECFEM2D with command: {cmd}")
            subprocess.run(cmd, shell=True, check=True,
                           stdout=subprocess.DEVNULL)

    def cleanup_xspecfem2d_run(self, choice=None):
        """
        Do some cleanup after running the SPECFEM2D binaries to make sure files
        are in the correct locations, and rename the OUTPUT_FILES directory so
        that it does not get overwritten by subsequent runs

        :type choice: str
        :param choice: Rename the OUTPUT_FILES directory with a suffix tag
            msut be 'INIT' or 'TRUE'. If None, will not rename but the
        """
        cd(self.workdir_paths.workdir)
        print("> Cleaning up after xspecfem2d, setting up for new run")

        # SPECFEM2D outputs its models in the DATA/ directory by default,
        # while SeisFlows expects this in the OUTPUT_FILES/ directory (which is
        # the default in SPECFEM3D)
        mv(glob.glob("DATA/*bin"), self.workdir_paths.output)

        if choice == "INIT":
            mv(self.workdir_paths.output, self.workdir_paths.model_init)
            # Create a new OUTPUT_FILES/ directory for TRUE run
            rm(self.workdir_paths.output)
            mkdir(self.workdir_paths.output)
        elif choice == "TRUE":
            mv(self.workdir_paths.output, self.workdir_paths.model_true)

    def setup_seisflows_working_directory(self):
        """
        Create and set the SeisFlows parameter file, making sure all required
        parameters are set correctly for this example problem
        """
        cd(self.cwd)

        self.sf.setup(force=True)  # Force will delete existing parameter file
        self.sf.configure()

        self.sf.par("ntask", self.ntask)  # default 3 sources for this example
        self.sf.par("materials", "elastic")  # how velocity model parameterized
        self.sf.par("density", "constant")  # update density or keep constant
        self.sf.par("format", "ascii")  # how to output synthetic seismograms
        self.sf.par("begin", 1)  # first iteration
        self.sf.par("end", self.niter)  # final iteration -- we will run 2
        self.sf.par("case", "synthetic")  # synthetic-synthetic inversion
        self.sf.par("attenuation", False)

        self.sf.par("specfem_bin", self.workdir_paths.bin)
        self.sf.par("specfem_data", self.workdir_paths.data)
        self.sf.par("model_init", self.workdir_paths.model_init)
        self.sf.par("model_true", self.workdir_paths.model_true)

    def finalize_specfem2d_par_file(self):
        """
        Last minute changes to get the SPECFEM2D Par_file in the correct format
        for running SeisFlows. Setting model type to read from .bin GLL files
        change number of stations etc.
        """
        cd(self.workdir_paths.data)
        self.sf.sempar("model", "gll")  # GLL so SPECFEM reads .bin files

    def run_sf_example(self):
        """
        Use subprocess to run the SeisFlows example we just set up
        """
        cd(self.cwd)
        subprocess.run("seisflows submit -f", check=False, shell=True)

    def main(self):
        """
        Setup the example and then optionally run the actual seisflows workflow
        """
        print(msg.cli("EXAMPLE SETUP", border="="))
        # Step 1: Download and configure SPECFEM2D, make binaries. Optional
        self.download_specfem2d()
        self.configure_specfem2d_and_make_binaries()
        # Step 2: Create a working directory and generate initial/final models
        self.create_specfem2d_working_directory()
        # Step 2a: Generate MODEL_INIT, rearrange consequent directory structure
        print(msg.cli("GENERATING INITIAL MODEL", border="="))
        self.setup_specfem2d_for_model_init()
        self.run_xspecfem2d_binaries()
        self.cleanup_xspecfem2d_run(choice="INIT")
        # Step 2b: Generate MODEL_INIT, rearrange consequent directory structure
        print(msg.cli("GENERATING TRUE/TARGET MODEL", border="="))
        self.setup_specfem2d_for_model_true()
        self.run_xspecfem2d_binaries()
        self.cleanup_xspecfem2d_run(choice="TRUE")
        # Step 3: Prepare Par_file and directory for MODEL_TRUE generation
        self.setup_seisflows_working_directory()
        self.finalize_specfem2d_par_file()
        print(msg.cli("COMPLETE EXAMPLE SETUP", border="="))
        # Step 4: Run the workflwo
        if self.run_example:
            print(msg.cli("RUNNING SEISFLOWS INVERSION WORKFLOW", border="="))
            self.run_sf_example()


if __name__ == "__main__":
    print(msg.ascii_logo_small)
    print(msg.cli(
        f"This is a [SPECFEM2D] [WORKSTATION] example, which will "
        f"run an inversion to assess misfit between two homogeneous halfspace "
        f"models with slightly different velocities. [3 events, 1 station, 2 "
        f"iterations]. The tasks involved include: ",
        items=["1. (optional) Download, configure, compile SPECFEM2D",
               "2. Set up a SPECFEM2D working directory",
               "3. Generate starting model from Tape2007 example",
               "4. Generate target model w/ perturbed starting model",
               "5. Set up a SeisFlows working directory",
               f"6. Run an inversion workflow"],
        header="seisflows example 1",
        border="=")
    )

    # Dynamically traverse sys.argv to get user-input command line. Cannot
    # use argparser here because we're being called by SeisFlows CLI tool which
    # is occupying argparser
    if len(sys.argv) > 1:
        sfex2d = SFExample2D()
        sfex2d.main()
