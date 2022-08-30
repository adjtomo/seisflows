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
from seisflows.tools import msg
from seisflows.tools.unix import cd
from seisflows.examples.sfexample2d import SFExample2D


class SFFwdEx2D(SFExample2D):
    """
    A class for running SeisFlows examples. Overloads Example 1 to take
    advantage of the default SPECFEM2D stuff, onyl changes the generation of
    MODEL TRUE, the number of stations, and the setup of the parameter file.
    """
    def __init__(self, ntask=None, nsta=None, method="run", specfem2d_repo=None,
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
        # Setting default values for ntask, niter, nsta here vvv
        super().__init__(ntask=ntask or 25, niter=1, nsta=nsta or 131,
                         method=method, specfem2d_repo=specfem2d_repo)

    def print_dialogue(self):
        """
        Print help/system dialogue message that explains the setup of th
        this workflow
        """
        print(msg.ascii_logo_small)
        print(msg.cli(
            f"This is a [SPECFEM2D] [WORKSTATION] example, which will "
            f"run forward simulations generate synthetic seismograms through "
            f"a given starting model. This example uses no preprocessing or "
            f"optimization modules"
            f"[{self.ntask} events, {self.nsta} stations, {self.niter} "
            f"iterations]. "
            f"The tasks involved include: ",
            items=["1. (optional) Download, configure, compile SPECFEM2D",
                   "2. Set up a SPECFEM2D working directory",
                   "3. Generate starting model from 'Tape2007' example",
                   "4. Set up a SeisFlows working directory",
                   "5. Run the forward simulation workflow"],
            header="seisflows example 2",
            border="=")
        )

    def setup_seisflows_working_directory(self):
        """
        Create and set the SeisFlows parameter file, making sure all required
        parameters are set correctly for this example problem
        """
        cd(self.cwd)

        print("> EX2: Setting SeisFlows parameters for Pyatao preprocessing")
        self.sf.setup(force=True)  # Force will delete existing parameter file
        self.sf.par("workflow", "forward")
        self.sf.par("preprocess", "null")
        self.sf.par("optimize", "null")
        self.sf.configure()

        self.sf.par("ntask", self.ntask)  # 3 sources for this example
        self.sf.par("materials", "elastic")  # how velocity model parameterized
        self.sf.par("density", False)  # update density or keep constant
        self.sf.par("data_format", "ascii")  # output synthetic seismograms
        self.sf.par("data_case", "synthetic")  # synthetic-synthetic inversion
        self.sf.par("attenuation", False)
        self.sf.par("components", "Y")

        self.sf.par("path_specfem_bin", self.workdir_paths.bin)
        self.sf.par("path_specfem_data", self.workdir_paths.data)
        self.sf.par("path_model_init", self.workdir_paths.model_init)

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
        # Step 3: Prepare Par_file and directory for MODEL_TRUE generation
        self.setup_seisflows_working_directory()
        self.finalize_specfem2d_par_file()
        print(msg.cli("COMPLETE EXAMPLE SETUP", border="="))
        # Step 4: Run the workflwo
        if self.run_example:
            print(msg.cli("RUNNING SEISFLOWS FORWARD WORKFLOW", border="="))
            self.run_sf_example()
            print(msg.cli("EXAMPLE COMPLETED SUCCESFULLY", border="="))
