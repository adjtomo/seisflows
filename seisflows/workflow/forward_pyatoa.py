#!/usr/bin/env python
"""
This is the subclass seisflows.workflow.forward_pyatoa

This is a main Seisflows class, it controls the main workflow.

This subclass inherits from the seisflows.workflow.inversion class
It includes additions to interact with the Python package Pyatoa, to perform
data collection, preprocessing, and misfit quantification steps
"""
import os
import sys
import time
import numpy as np
from glob import glob

from seisflows.tools.err import ParameterError
from seisflows.tools.tools import exists
from seisflows.config import custom_import

from pyatoa import Pyaflowa

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']
optimize = sys.modules['seisflows_optimize']
preprocess = sys.modules['seisflows_preprocess']


class ForwardPyatoa(custom_import("workflow", "forward")):
    """
    Forward solver subclass, additional support for Pyatoa integration
    """

    def check(self):
        """
        Checks parameters and paths
        """
        # Run Base class checks
        super(ForwardPyatoa, self).check()

        # Signifiy if data-synth. or synth.-synth. case
        if "CASE" not in PAR:
            raise ParameterError(PAR, "CASE")

        if "RESUME_FROM" not in PAR:
            setattr(PAR, "RESUME_FROM", None)

        # Synthetic-synthetic examples require a true model to create the "data"
        if PAR.CASE == "Synthetic" and not exists(PATH.MODEL_TRUE):
            raise Exception()

        # Pyatoa specific paths
        if "PYATOA_IO" not in PATH:
            raise ParameterError(PATH, "PYATOA_IO")

    def setup(self):
        """
        Overwrite seisflows.workflow.inversion.setup()

        Lays groundwork for inversion, removes need for a preprocess.setup(),
        instead runs a setup for Pyatoa
        """
        print("\tInitializing Pyaflowa")
        self.pyaflowa = Pyaflowa(par=vars(PAR), paths=vars(PATH))

        print("\tPreparing initial model", end="... ")
        self.stopwatch("set")
        system.run("solver", "setup")
        self.stopwatch("time")

        self.prepare_model()

    def initialize(self):
        """
        Overwrite seisflows.workflow.inversion.intialize()

        Generates synthetics via a forward simulation, calculates misfits
        for the forward simulation. Writes misfit for use in optimization.

        Breaks apart the function evaluatoin into a forward simulation and a
        function evaluation
        """
        print("INITIALIZE")
        suffix_ = "new"
        path_ = PATH.SCRATCH

        print("\tRunning forward simulation", end="... ")
        self.stopwatch("set")
        system.run("solver", "eval_fwd", path=path_)
        self.stopwatch("time")

        print("\tQuantifying misfit", end="... ")
        self.stopwatch("set")
        self.pyaflowa.set(iteration=0, step=0)
        system.run_ancil("solver", "eval_func", pyaflowa=self.pyaflowa)
        self.stopwatch("time")

    def write_misfit(self, suffix):
        """
        Overwrites seisflows.workflow.inversion.write_misfit()

        Total misfit defined by Tape et al 2007: total_misfit=misfit/number_srcs

        Writes misfit in format expected by nonlinear optimization library
        Waits for all instances of Pyatoa to finish writing their misfit

        :type suffix: str
        :param suffix: suffix to write the misfit with, e.g. `f_new`
        """
        src = os.path.join(PATH.PYATOA_IO, "data", "misfits", "*")
        dst = f"f_{suffix}"

        # Make sure pyaflowa has written all misfits before proceeding
        # Wait as long as it would reasonably take for a task to finish
        waited = 0
        wait_time_s = 10
        while len(glob(src)) < PAR.NTASK:
            time.sleep(wait_time_s)
            waited += wait_time_s
            if waited >= (PAR.ANCIL_TASKTIME * 60) / 2:
                print("\t\tworkflow.inversion_pyatoa.write_misfit: waited too "
                      "long for misfits to be written. exiting.")
                sys.exit(-1)

        # Sum up individual event misfits, remove them afterwards
        total_misfit = 0
        for event in glob(src):
            total_misfit += np.loadtxt(event)
            os.remove(event)

        # Save the total misfit
        total_misfit = total_misfit / PAR.NTASK
        optimize.savetxt(dst, total_misfit)




