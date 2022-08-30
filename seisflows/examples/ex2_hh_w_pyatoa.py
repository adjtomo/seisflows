#!/usr/bin/env python3
"""
                SEISFLOWS SPECFEM2D WORKSTATION EXAMPLE 2

This example will run N iterations of an inversion to assess misfit between
a homogeneous halfspace model and a checkerboard model using X events and
Y receivers. N, X and Y are all user-selectable.

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
    def __init__(self, ntask=None, niter=None, nsta=None, method="run",
                 specfem2d_repo=None):
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
        # Setting default values for ntask, niter, nsta here vvv
        super().__init__(ntask=ntask or 4, niter=niter or 2, nsta=nsta or 32,
                         method=method, specfem2d_repo=specfem2d_repo)

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

    def print_dialogue(self):
        """
        Print help/system dialogue message that explains the setup of th
        this workflow
        """
        print(msg.ascii_logo_small)
        print(msg.cli(
            f"This is a [SPECFEM2D] [WORKSTATION] example, which will "
            f"run an inversion to assess misfit between a starting homogeneous "
            f"halfspace model and a target checkerboard model. This "
            f"example problem uses the [PYAFLOWA] preprocessing "
            f"module and the [LBFGS] optimization algorithm. "
            f"[{self.ntask} events, {self.nsta} stations, {self.niter} "
            f"iterations]. "
            f"The tasks involved include: ",
            items=["1. (optional) Download, configure, compile SPECFEM2D",
                   "2. Set up a SPECFEM2D working directory",
                   "3. Generate starting model from 'Tape2007' example",
                   "4. Generate target model w/ perturbed starting model",
                   "5. Set up a SeisFlows working directory",
                   "6. Run the inversion workflow"],
            header="seisflows example 2",
            border="=")
        )

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


