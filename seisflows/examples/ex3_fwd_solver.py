#!/usr/bin/env python3
"""
                SEISFLOWS SPECFEM2D WORKSTATION EXAMPLE 3

This example will run a number of forward simulations and misfit quantification.
This is useful for generating a large number of synthetics through a given model

.. note::
    See Example 1 docstring for more information

.. rubric::
    $ seisflows examples run 3
"""
import os
from seisflows.tools import msg
from seisflows.tools.unix import cd, ln, rm
from seisflows.examples.sfexample2d import SFExample2D


class SFFwdEx2D(SFExample2D):
    """
    A class for running SeisFlows examples. Overloads Example 1 to take
    advantage of the default SPECFEM2D stuff, only removes the generation of
    MODEL TRUE and makes the parameter file more bare bones.
    """
    def __init__(self, ntask=None, nsta=None, nproc=None, method="run",
                 specfem2d_repo=None, with_mpi=False, mpiexec="mpirun",
                 **kwargs):
        """
        Overloads init of the base problem

        :type ntask: int
        :param ntask: number of events to use in inversion, between 1 and 25.
            defaults to 3
        :type nsta: int
        :param nsta: number of stations to include in inversion, between 1 and
            131
        :type method: str
        :param method: method for running the example problem, can be:
            * 'run': setup and run the example problem
            * 'setup': only setup the example problem, do not `submit` job
        :type specfem2d_repo: str
        :param specfem2d_repo: path to the SPECFEM2D directory which should
            contain binary executables. If not given, SPECFEM2D will be
            downloaded configured and compiled automatically.
        """
        # We set defaults here because `seisflows examples` may input these
        # values as NoneType which would override __init__ defaults.
        super().__init__(ntask=ntask or 10, nsta=nsta or 25, nproc=nproc or 1,
                         niter=1, method=method, specfem2d_repo=specfem2d_repo,
                         with_mpi=with_mpi, mpiexec=mpiexec)

        self._modules = {
            "workflow": "forward",
            "preprocess": "null",
            "optimize": "null"
        }

        self._parameters["export_traces"] = True
        self._parameters["path_model_true"] = "null"  # overload default par.

    def print_dialogue(self):
        """
        Print help/system dialogue message that explains the setup of th
        this workflow
        """
        print(msg.ascii_logo_small)
        print(msg.cli(
            f"This is a [SPECFEM2D] [WORKSTATION] example, which will run "
            f"forward simulations to generate synthetic seismograms through "
            f"a homogeneous halfspace starting model. This example uses no "
            f"preprocessing or optimization modules. "
            f"[{self.ntask} events, {self.nsta} stations] "
            f"The tasks involved include: ",
            items=["1. (optional) Download, configure, compile SPECFEM2D",
                   "2. [Setup] a SPECFEM2D working directory",
                   "3. [Setup] starting model from 'Tape2007' example",
                   "4. [Setup] a SeisFlows working directory",
                   "5. [Run] the forward simulation workflow"],
            header="seisflows example 3",
            border="=")
        )

    def setup_specfem2d_for_model_true(self):
        """
        Overwrites MODEL TRUE creation from EX1. The same as in Example 2

        Make some adjustments to the parameter file to create the final velocity
        model. This function assumes it is running from inside the
        SPECFEM2D/DATA directory
        """
        cd(self.workdir_paths.data)
        assert(os.path.exists("Par_file")), f"I cannot find the Par_file!"

        print("> Updating SPECFEM2D to set checkerboard model as current model")
        self.sf.sempar("model", "legacy")  # read model_velocity.dat_checker
        rm("proc000000_model_velocity.dat_input")
        for i in range(self.nproc):
            ln("model_velocity.dat_checker",
               f"proc00000{i}_model_velocity.dat_input")


    def main(self):
        """
        Setup the example and then optionally run the actual seisflows workflow
        Mostly the same as Example 1 main() except it does not generate
        MODEL_TRUE, and instead sets MODEL_TRUE as the starting model.
        """
        print(msg.cli("EXAMPLE SETUP", border="="))

        # Step 1: Download and configure SPECFEM2D, make binaries. Optional
        self.download_specfem2d()
        self.configure_specfem2d()
        self.make_specfem2d_executables()
        # Step 2: Create a working directory and generate initial/final models
        self.create_specfem2d_working_directory()
        # Step 2a: Generate MODEL_INIT, rearrange consequent directory structure
        print(msg.cli("GENERATING INITIAL MODEL", border="="))
        self.setup_specfem2d_for_model_init()  # setup SPECFEM run directory
        self.setup_specfem2d_for_model_true()  # Use checkerboard model
        self.run_xspecfem2d_binaries()
        self.cleanup_xspecfem2d_run(choice="INIT")
        # Step 3: Prepare Par_file and directory for MODEL_TRUE generation
        self.setup_seisflows_working_directory()
        self.finalize_specfem2d_par_file()
        print(msg.cli("COMPLETE EXAMPLE SETUP", border="="))
        # Step 4: Run the workflwo
        if self.run_example:
            self.run_sf_example()
