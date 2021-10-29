#!/usr/bin/env python
"""
This is a subclass seisflows.workflow.InversionMaui
"""
import sys

from seisflows3.config import custom_import
from seisflows3.tools.err import ParameterError

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]

system = sys.modules["seisflows_system"]
solver = sys.modules["seisflows_solver"]
optimize = sys.modules["seisflows_optimize"]
preprocess = sys.modules["seisflows_preprocess"]
postprocess = sys.modules["seisflows_postprocess"]


class ThriftyMaui(custom_import("workflow", "thrifty_inversion")):
    """
    Waveform thrify inversion class specifically for running jobs on the
    New Zealand HPC cluster Maui.

    On Maui, Anaconda is only available on an ancillary cluster, Maui_ancil,
    so jobs involving the preprocessing module must be called through a
    separate system run call.
    """
    def check(self):
        """
        Ensure that the correct submodules are specified, otherwise
        this workflow won't function properly.
        """
        super().check()

        if "MAUI" not in PAR.SYSTEM.upper():
            raise ParameterError()

        if "MAUI" not in PAR.SOLVER.upper():
            raise ParameterError()

    def setup(self):
        """
        Lays groundwork for inversion by running setup() functions for the
        involved sub-modules, and generating synthetic true data if necessary,
        and generating the pre-requisite database files. Should only be run once
        at the iteration 1
        """
        # Set up all the requisite modules
        print("SETUP")
        preprocess.setup()
        postprocess.setup()
        optimize.setup()

        # Run the setup in serial to reduce unnecessary job submissions
        # Needs to be split up into multiple system calls
        solver.initialize_solver_directories()

        if PAR.CASE.upper() == "SYNTHETIC":
            system.run_single("solver", "setup", model="true")
            system.run("solver", "generate_data")

        system.run_single("solver", "setup", model="init")

    def evaluate_function(self, path, suffix):
        """
        Performs forward simulation, and evaluates the objective function.

        Differs from Inversion.evaluate_function() as it splits the forward
        problem and misfit quantification into two separate system calls,
        rather than a single system call.

        :type path: str
        :param path: path in the scratch directory to use for I/O
        :type suffix: str
        :param suffix: suffix to use for I/O
        """
        print("EVALUATE FUNCTION\n\tRunning forward simulation")
        self.write_model(path=path, suffix=suffix)
        system.run("solver", "eval_fwd", path=path)
        print("\tEvaluating misfit")
        system.run_ancil("solver", "eval_misfit", path=path)
        self.write_misfit(path=path, suffix=suffix)

