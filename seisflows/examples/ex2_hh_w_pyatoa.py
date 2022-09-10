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
    def __init__(self, ntask=None, niter=None, nsta=None, nproc=None,
                 event_id=None, method="run", specfem2d_repo=None,
                 with_mpi=False, mpiexec="mpirun", **kwargs):
        """
        Overload init and attempt to import Pyatoa before running example.

        :type ntask: int
        :param ntask: number of events to use in inversion, between 1 and 25.
            defaults to 3
        :type niter: int
        :param niter: number of iterations to run. defaults to 2
        :type nsta: int
        :param nsta: number of stations to include in inversion, between 1 and
            131
        :type event_id: str
        :param event_id: allow user to choose a specific event ID from the
            example problem. Must match source files in SPECFEM2D example DATA/
            directory. Overwrites `ntask` to be 1
        :type specfem2d_repo: str
        :param specfem2d_repo: path to the SPECFEM2D directory which should
            contain binary executables. If not given, SPECFEM2D will be
            downloaded configured and compiled automatically.
        """
        # Because of how the checkerboard model is defined (as a .dat file),
        # we cannot run with nproc > 1 because the checkerboard model will be
        # partititioned incorrectly. MPI with nproc==1 still okay.
        if nproc is not None:
            assert(nproc == 1), \
                f"Example #2 can only be run as a serial (nproc==1) task"

        # We set defaults here because `seisflows examples` may input these
        # values as NoneType which would override __init__ defaults.
        super().__init__(ntask=ntask or 4, niter=niter or 2, nsta=nsta or 32,
                         event_id=event_id, method=method, 
                         specfem2d_repo=specfem2d_repo, with_mpi=with_mpi,
                         mpiexec=mpiexec)

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

        # Define the main SeisFlows modules, which will be different to the
        # base example
        self._modules = {
            "workflow": "inversion",
            "preprocess": "pyaflowa",
            "optimize": "LBFGS",
        }

        # Adjust the existing parameter list
        self._parameters["smooth_h"] = 5000.
        self._parameters["smooth_v"] = 5000.
        self._parameters["pyflex_preset"] = "null"  # no windowing in Pyaflowa

        # Pyaflowa preprocessing parameters
        self._parameters["unit_output"] = "DISP"
        self._parameters["min_period"] = 10.  # filter bounds define windows
        self._parameters["max_period"] = 200.

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
            f"module and the [LBFGS] optimization algorithm. In this example, "
            f"windowing in Pyaflowa is turnd OFF."
            f"[{self.ntask} events, {self.nsta} stations, {self.niter} "
            f"iterations]. "
            f"The tasks involved include: ",
            items=["1. (optional) Download, configure, compile SPECFEM2D",
                   "2. [Setup] a SPECFEM2D working directory",
                   "3. [Setup] starting model from 'Tape2007' example",
                   "4. [Setup] target model w/ perturbed starting model",
                   "5. [Setup] a SeisFlows working directory",
                   "6. [Run] the inversion workflow"],
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

        print("> Updating SPECFEM2D to set checkerboard model as current model")
        self.sf.sempar("model", "legacy")  # read model_velocity.dat_checker
        rm("proc000000_model_velocity.dat_input")
        ln("model_velocity.dat_checker", "proc000000_model_velocity.dat_input")

