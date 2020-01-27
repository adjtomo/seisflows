#!/usr/bin/env python
"""
This is the subclass seisflows.workflow.InversionPyatoa

This is a main Seisflows class, it controls the main workflow.

This subclass inherits from the seisflows.workflow.inversion class
It includes additions to interact with the Python package Pyatoa, to perform
data collection, preprocessing, and misfit quantification steps
"""
import os
import sys
import time
import glob
import numpy as np

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
postprocess = sys.modules['seisflows_postprocess']


class InversionPyatoa(custom_import('workflow', 'inversion')):
    """
    Waveform inversion subclass, additional support for Pyatoa integration
    
    Overwrites following functions from seisflows.workflow.inversion.Inversion
        
        check, setup, initialize, evaluate_function, finalize, write_misfit
    """
    def check(self):
        """
        Checks parameters and paths
        """
        # Run Base class checks
        super(InversionPyatoa, self).check()

        # Signifiy if data-synth. or synth.-synth. case
        if 'CASE' not in PAR:
            raise ParameterError(PAR, 'CASE')

        # Synthetic-synthetic examples require a true model to create the 'data'
        if PAR.CASE == 'Synthetic' and not exists(PATH.MODEL_TRUE):
            raise Exception()

        # Pyatoa specific paths
        if 'PYATOA_IO' not in PATH:
            raise ParameterError(PATH, 'PYATOA_IO')

    def setup(self):
        """
        Overwrite seisflows.workflow.inversion.setup()

        Lays groundwork for inversion, removes need for a preprocess.setup(),
        instead runs a setup for Pyatoa
        """
        if optimize.iter == 1:
            print("SETUP\n\tPerforming module setup")
            postprocess.setup()
            optimize.setup()

            print("\tInitializing Pyaflowa")
            self.pyaflowa = Pyaflowa(par=vars(PAR), paths=vars(PATH))

            print("\tPreparing initial model", end="... ")
            self.stopwatch("set")
            system.run("solver", "setup")
            self.stopwatch("time")

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
        path_ = PATH.GRAD

        self.write_model(path=path_, suffix=suffix_)

        print("\tRunning forward simulation", end="... ")
        self.stopwatch("set")
        system.run("solver", "eval_fwd", path=path_)
        self.stopwatch("time")

        print("\tQuantifying misfit", end="... ")
        self.stopwatch("set")
        self.pyaflowa.set(iteration=optimize.iter, step=0, misfits={})
        system.run_ancil("solver", "eval_func", pyaflowa=self.pyaflowa)
        self.stopwatch("time")

        print("\tWriting misfit")
        self.write_misfit(suffix=suffix_)

    def evaluate_function(self):
        """
        Overwrite seisflows.workflow.inversion.evaluate_function()

        This is the same function as intialize, however the input and output
        paths are different, to signify that this is executed at a different
        section of the workflow
        """
        suffix_ = "try"
        path_ = PATH.FUNC

        self.write_model(path=path_, suffix=suffix_)

        print("\tRunning forward simulation", end="... ")
        self.stopwatch("set")
        system.run("solver", "eval_fwd", path=path_)
        self.stopwatch("time")

        print("\tQuantifying misfit", end="... ")
        self.stopwatch("set")
        self.pyaflowa.set(iteration=optimize.iter,
                          step=optimize.line_search.step_count + 1,
                          misfits={}
                          )
        system.run_ancil("solver", "eval_func", pyaflowa=self.pyaflowa)
        self.stopwatch("time")

        self.write_misfit(suffix=suffix_)

    def finalize(self):
        """
        Overwrites seisflows.workflow.inversion.finalize()

        Adds finalization of Pyaflowa and saves results from current iteration
        """
        super(InversionPyatoa, self).finalize()

        # Finalize Pyatoa for the given iteration and step count
        print("\tFinalizing Pyatoa", end="... ")
        self.stopwatch("set")
        self.pyaflowa.set(iteration=optimize.iter,
                          step=optimize.line_search.step_count)
        self.pyaflowa.finalize()
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
        dst = f"f_{suffix}"

        # Make sure pyaflowa has written all misfits before proceeding
        waited = 0
        while len(self.pyaflowa.misfits) < PAR["NSRC"]:
            time.sleep(10)
            waited += 1
            if waited >= 10:
                print("seisflows.workflow.inversion_pyatoa waited too long for "
                      "misfits to be written, exiting.")
                sys.exit(-1)

        # Save the total misfit
        total_misfit = sum(self.pyaflowa.misfits.values()) / PAR["NSRC"]
        optimize.savetxt(dst, total_misfit)




