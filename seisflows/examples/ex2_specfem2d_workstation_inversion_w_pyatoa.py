#!/usr/bin/env python3
"""
                SEISFLOWS SPECFEM2D WORKSTATION EXAMPLE 2

This example will run two iterations of an inversion to assess misfit between
a homogeneous halfspace model and a checkerboard model using 2 events and
5 receivers.

.. note::
    See Example 1 docstring for more information

.. rubric::
    $ seisflows examples run 2
"""
import os
import sys

from seisflows.tools import msg
from seisflows.tools.unix import cd, rm, ln
from seisflows.examples.sfexample2d import SFExample2D


class SFPyatoaEx2D(SFExample2D):
    """
    A class for running SeisFlows examples. Overloads Example 1 to take
    advantage of the default SPECFEM2D stuff, onyl changes the generation of
    MODEL TRUE, the number of stations, and the setup of the parameter file.
    """
    def __init__(self, ntask=2, niter=2, nsta=5, specfem2d_repo=None):
        """
        Overload init and attempt to import Pyatoa before running example,
        overload the default number of tasks to 2, and add a new init parameter
        `nsta` which chooses the number of stations, between 1 and 132

        :type ntask: int
        :param ntask: number of events to use in inversion, between 1 and 25.
            defaults to 3
        :type niter: int
        :param niter: number of iterations to run. defaults to 2
        :type nsta: int
        :param nsta: number of stations to include in inversion, between 1 and
            131
        :type specfem2d_repo: str
        :param specfem2d_repo: path to the SPECFEM2D directory which should
            contain binary executables. If not given, SPECFEM2D will be
            downloaded configured and compiled automatically.
        """
        super().__init__(ntask=ntask, niter=niter,
                         specfem2d_repo=specfem2d_repo)
        self.nsta = nsta
        # -1 because it represents index but we need to talk in terms of count
        assert(1 <= self.nsta <= 131), \
            f"number of stations must be between 1 and 131, not {self.nsta}"
        # Make sure that Pyatoa has been installed before running
        try:
            import pyatoa
        except ModuleNotFoundError:
            print(msg.cli("Module Pyatoa not found but is required for this "
                          "example. Please install Pyatoa and rerun this "
                          "example.", header="module not found error",
                          border="=")
                  )
            sys.exit(-1)

    def setup_specfem2d_for_model_true(self):
        """
        Overwrites MODEL TRUE creation from EX1

        Make some adjustments to the parameter file to create the final velocity
        model. This function assumes it is running from inside the
        SPECFEM2D/DATA directory
        """
        cd(self.workdir_paths.data)
        assert(os.path.exists("Par_file")), f"I cannot find the Par_file!"

        print("> EX: Updating SPECFEM2D to set checkerboard model as "
              "MODEL_TRUE")
        self.sf.sempar("model", "legacy")  # read model_velocity.dat_checker
        rm("proc000000_model_velocity.dat_input")
        ln("model_velocity.dat_checker", "proc000000_model_velocity.dat_input")

    def setup_seisflows_working_directory(self):
        """
        Create and set the SeisFlows parameter file, making sure all required
        parameters are set correctly for this example problem
        """
        cd(self.cwd)

        print("> EX2: Setting SeisFlows parameters for Pyatao preprocessing")
        self.sf.setup(force=True)  # Force will delete existing parameter file
        self.sf.par("workflow", "inversion")
        self.sf.par("preprocess", "pyaflowa")
        self.sf.par("optimize", "LBFGS")
        self.sf.configure()

        self.sf.par("end", 1)  # only 1 iteration
        self.sf.par("ntask", self.ntask)  # 3 sources for this example
        self.sf.par("materials", "elastic")  # how velocity model parameterized
        self.sf.par("density", False)  # update density or keep constant
        self.sf.par("data_format", "ascii")  # output synthetic seismograms
        self.sf.par("data_case", "synthetic")  # synthetic-synthetic inversion
        self.sf.par("attenuation", False)
        self.sf.par("components", "Y")

        # PYATOA preprocessing parameters
        self.sf.par("unit_output", "DISP")
        self.sf.par("min_period", 10)  # filter bounds define window selection
        self.sf.par("max_period", 200)
        # self.sf.par("pyflex_preset", "")  # To turn off windowing completely

        self.sf.par("path_specfem_bin", self.workdir_paths.bin)
        self.sf.par("path_specfem_data", self.workdir_paths.data)
        self.sf.par("path_model_init", self.workdir_paths.model_init)
        self.sf.par("path_model_true", self.workdir_paths.model_true)

    def finalize_specfem2d_par_file(self):
        """
        Final changes to the SPECFEM2D Par_file before running SeisFlows.
        Par_file will be used to control all the child specfem2d directories.
        Need to tell them to read models from .bin files, and to use existing
        station files rather than create them from the Par_file
        """
        print("> EX2: Finalizing SPECFEM2D Par_file for SeisFlows inversion")

        cd(self.workdir_paths.data)
        self.sf.sempar("model", "gll")  # GLL so SPECFEM reads .bin files
        self.sf.sempar("use_existing_stations", ".true.")  # Use STATIONS file
        # Assign STATIONS_checker file which has 132 stations
        rm("STATIONS")

        # Only write the first 10 lines to get 10 stations in inversion
        with open("STATIONS_checker", "r") as f:
            lines = f.readlines()

        print(f"> EX2: Using {self.nsta} stations in this inversion workflow")
        with open("STATIONS", "w") as f:
            f.writelines(lines[:self.nsta])


if __name__ == "__main__":
    print(msg.ascii_logo_small)
    print(msg.cli(
        f"This is a [SPECFEM2D] [WORKSTATION] example, which will "
        f"run an inversion to assess misfit between a homogeneous halfspace  "
        f"and checkerboard model using Pyatoa for misfit quantification "
        f"[2 events, 5 stations, 1 iterations]. The tasks involved include: ",
        items=["1. (optional) Download, configure, compile SPECFEM2D",
               "2. Set up a SPECFEM2D working directory",
               "3. Generate starting model from Tape2007 example",
               "4. Generate target model w/ perturbed starting model",
               "5. Set up a SeisFlows working directory",
               f"6. Run an inversion workflow. The line search is expected to "
               f"attempt 4 evaluations (i01s04)"],
        header="seisflows example 2",
        border="=")
    )

    # Dynamically traverse sys.argv to get user-input command line. Cannot
    # use argparser here because we're being called by SeisFlows CLI tool which
    # is occupying argparser
    if len(sys.argv) > 1:
        _, _, specfem2d_repo = sys.argv
        sfex2d = SFPyatoaEx2D(specfem2d_repo=specfem2d_repo)
        sfex2d.main()