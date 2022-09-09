#!/usr/bin/env python3
"""
                SEISFLOWS SPECFEM2D WORKSTATION EXAMPLE 4

This example mimics Example 3 (mass forward simulations) but uses the
parallelized version of SPECFEM2D to test out multi-core functionality.

.. note::
    See Example 1 docstring for more information

.. rubric::
    $ seisflows examples run 3
"""
import os
import subprocess
from seisflows.tools import msg
from seisflows.tools.unix import cd, rm, ln
from seisflows.examples.sfexample2d import SFExample2D


class SFMultiCoreEx2D(SFExample2D):
    """
    A class for running SeisFlows examples. Overloads Example 1 to take
    advantage of the default SPECFEM2D stuff, onyl changes the generation of
    MODEL TRUE, the number of stations, and the setup of the parameter file.
    """
    def __init__(self, ntask=None, nsta=None, nproc=None,
                 method="run", specfem2d_repo=None, **kwargs):
        """
        Overloads init of the base problem

        :type ntask: int
        :param ntask: number of events to use in inversion, between 1 and 25.
            defaults to 3
        :type nsta: int
        :param nsta: number of stations to include in inversion, between 1 and
            131
        :type nproc: int
        :param nproc: number of processors to be sent to MPI executable
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
        super().__init__(ntask=ntask or 10, nsta=nsta or 25, nproc=nproc or 4,
                         niter=1, method=method, specfem2d_repo=specfem2d_repo)

        self.mpiexec = "mpirun"
        self._parameters["nproc"] = self.nproc
        self._parameters["mpiexec"] = self.mpiexec

        self._modules = {
            "workflow": "forward",
            "preprocess": "null",
            "optimize": "null",
        }

        self._parameters["export_traces"] = True
        self._parameters["path_model_true"] = "null"  # overload default par.

        # Overwrite configure cmd to get MPI
        self._configure_cmd = \
            "./configure FC=gfortran CC=gcc MPIF90=mpif90 --with-mpi"


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
            f"preprocessing or optimization modules. This example uses MPI to "
            f"run the external solver, SPECFEM2D."
            f"[{self.ntask} events, {self.nsta} stations] "
            f"The tasks involved include: ",
            items=["1. (optional) Download, configure, compile "
                   "SPECFEM2D w/ MPI",
                   "2. [Setup] a SPECFEM2D working directory",
                   "3. [Setup] starting model from 'Tape2007' example",
                   "4. [Setup] a SeisFlows working directory",
                   "5. [Run] the forward simulation workflow"],
            header="seisflows example 4",
            border="=")
        )

    def run_xspecfem2d_binaries(self):
        """
        Runs the xmeshfem2d and then xspecfem2d binaries using subprocess and then
        do some cleanup to get files in the correct locations. Assumes that we
        can run the binaries directly with './'
        """
        cd(self.workdir_paths.workdir)

        mpicmd = f"{self.mpiexec} -n {self.nproc}"
        cmd_mesh = f"{mpicmd} bin/xmeshfem2D > OUTPUT_FILES/mesher.log.txt"
        cmd_spec = f"{mpicmd} bin/xspecfem2D > OUTPUT_FILES/solver.log.txt"

        for cmd in [cmd_mesh, cmd_spec]:
            print(f"Running SPECFEM2D with command: {cmd}")
            subprocess.run(cmd, shell=True, check=True,
                           stdout=subprocess.DEVNULL)

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
        self.run_xspecfem2d_binaries()
        self.cleanup_xspecfem2d_run(choice="INIT")
        # Step 3: Prepare Par_file and directory for MODEL_TRUE generation
        self.setup_seisflows_working_directory()
        self.finalize_specfem2d_par_file()
        print(msg.cli("COMPLETE EXAMPLE SETUP", border="="))
        # Step 4: Run the workflwo
        if self.run_example:
            self.run_sf_example()
