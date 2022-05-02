#!/usr/bin/env python3
"""
                SEISFLOWS3 SPECFEM2D WORKSTATION EXAMPLE 2

This example will run two iterations of an inversion to assess misfit between
a homogeneous halfspace model and a checkerboard model using 3 events and
132 receivers.

.. note::
    See Example 1 docstring for more information

.. rubric::
    $ seisflows examples run 2
"""
import os
import sys
import glob
import shutil
import subprocess
import numpy as np

from seisflows3.tools import msg
from seisflows3.config import Dict
from seisflows3.seisflows import SeisFlows
from seisflows3.tools.unix import cd, cp, rm, ln, mv, mkdir
from seisflows3.scripts.examples.ex1_specfem2d_workstation_inversion import \
    SF3Example2D


class SF3PyatoaEx2D(SF3Example2D):
    """
    A class for running SeisFlows3 examples. Overloads Example 1 to take
    advantage of the default SPECFEM2D stuff, onyl changes the generation of
    MODEL TRUE, the number of stations, and the setup of the parameter file.
    """
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
        Create and set the SeisFlows3 parameter file, making sure all required
        parameters are set correctly for this example problem
        """
        cd(self.cwd)

        print("> EX2: Setting SeisFlows3 parameters for Pyatao preprocessing")
        self.sf.setup(force=True)  # Force will delete existing parameter file
        self.sf.par("preprocess", "pyatoa")
        self.sf.configure()

        self.sf.par("end", 1)  # only 1 iteration
        self.sf.par("ntask", self.ntask)  # we will be using 3 sources for this example
        self.sf.par("materials", "elastic")  # how the velocity model is parameterized
        self.sf.par("density", "constant")  # update density or keep constant
        self.sf.par("nt", 5000)  # set by SPECFEM2D Par_file
        self.sf.par("dt", .06)  # set by SPECFEM2D Par_file
        self.sf.par("f0", 0.084)  # set by SOURCE file
        self.sf.par("format", "ascii")  # how to output synthetic seismograms
        self.sf.par("case", "synthetic")  # synthetic-synthetic inversion
        self.sf.par("attenuation", False)
        self.sf.par("components", "Y")

        # PYATOA preprocessing parameters
        self.sf.par("unit_output", "DISP")
        self.sf.par("min_period", 10)  # filter bounds define window selection
        self.sf.par("max_period", 200)
        self.sf.par("start_pad", 48)  # T0 set in Par_file
        self.sf.par("end_pad", 5000 * .06)  # nt * dt defined by Par_file
        # self.sf.par("pyflex_preset", "")  # To turn off windowing completely

        self.sf.par("specfem_bin", self.workdir_paths.bin)
        self.sf.par("specfem_data", self.workdir_paths.data)
        self.sf.par("model_init", self.workdir_paths.model_init)
        self.sf.par("model_true", self.workdir_paths.model_true)

    def finalize_specfem2d_par_file(self):
        """
        Final changes to the SPECFEM2D Par_file before running SeisFlows. Par_file
        will be used to control all the child specfem2d directories. Need to tell
        them to read models from .bin files, and to use existing station files
        rather than create them from the Par_file
        rather than create them from the Par_file
        """
        print("> EX2: Finalizing SPECFEM2D Par_file for SeisFlows3 inversion")

        cd(self.workdir_paths.data)
        self.sf.sempar("model", "gll")  # GLL so SPECFEM reads .bin files
        self.sf.sempar("use_existing_stations", ".true.")  # Use STATIONS file
        # Assign STATIONS_checker file which has 132 stations
        rm("STATIONS")
        ln("STATIONS_checker", "STATIONS")


if __name__ == "__main__":
    print(msg.ascii_logo_small)
    print(msg.cli(
        f"This is a [SPECFEM2D] [WORKSTATION] example, which will "
        f"run an inversion to assess misfit between a homogeneous halfspace  "
        f"and checkerboard model using Pyatoa for misfit quantification "
        f"[3 events, 132 station, 1 iterations]. The tasks involved include: ",
        items=["1. (optional) Download, configure, compile SPECFEM2D",
               "2. Set up a SPECFEM2D working directory",
               "3. Generate starting model from Tape2007 example",
               "4. Generate target model w/ perturbed starting model",
               "5. Set up a SeisFlows3 working directory",
               f"6. Run an inversion workflow"],
        header="seisflows3 example 2",
        border="=")
    )

    # Dynamically traverse sys.argv to get user-input command line. Cannot
    # use argparser here because we're being called by SeisFlows CLI tool which
    # is occupying argparser
    if len(sys.argv) > 1:
        sf3ex2d = SF3PyatoaEx2D()
        sf3ex2d.main()
