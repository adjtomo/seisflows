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

    def main(self):
        """
        !!! This function controls the main workflow !!!

        Overwrites seisflows.workflow.inversion.main.
        Carries out seismic inversion
        """
        # Make the workflow a list of functions that can be called dynamically
        flow = [self.initialize,
                self.evaluate_gradient,
                self.write_gradient,
                self.compute_direction,
                self.line_search,
                self.finalize,
                self.clean
                ]

        self.stopwatch("set")
        print(f"BEGINNING WORKFLOW AT {self.stopwatch()}")
        optimize.iter = PAR.BEGIN
        print(f"{optimize.iter} <= {PAR.END}")

        # Allow workflow resume from a given mid-workflow location
        if PAR.RESUME_FROM:
            self.resume_from(flow)
        else:
            # First-time intialization of the workflow
            self.setup()

        # Run the workflow until PAR.END
        while optimize.iter <= PAR.END:
            print(f"ITERATION {optimize.iter}")
            for func in flow:
                func()
            print(f"FINISHED ITERATION {optimize.iter} AT {self.stopwatch()}\n")
            self.stopwatch("time")
            optimize.iter += 1

    def resume_from(self, flow):
        """
        Resume the workflow from a given function. Unique function
    
        :type flow: list
        :param flow: list of functions which comprise the full workflow
        """ 
        # Determine the index that corresponds to the resume function named
        for i, func in enumerate(flow):
            if func.__name__ == PAR.RESUME_FROM:
                resume_idx = i
                break
        else:
            print("PAR.RESUME_FROM does not correspond to any workflow "
                  "functions. Exiting...")
            sys.exit(-1)
        print(f"RESUME ITERATION {optimize.iter} (from function "
              f"{flow[resume_idx].__name__})")
        
        for func in flow[resume_idx:]:
            func()
        print(f"FINISHED ITERATION {optimize.iter} AT {self.stopwatch()}\n")
        optimize.iter += 1

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
        self.pyaflowa.set(iteration=optimize.iter, step=0)
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
                          step=optimize.line_search.step_count + 1)
        system.run_ancil("solver", "eval_func", pyaflowa=self.pyaflowa)
        self.stopwatch("time")

        self.write_misfit(suffix=suffix_)

    def evaluate_gradient(self):
        """
        Overwrites seisflows.workflow.inversion.evaluate_gradient()        

        Performs adjoint simulation to evaluate gradient of the misfit function
        """
        print("EVALUATE GRADIENT\n\tRunning adjoint simulation", end="... ")

        self.stopwatch("set")
        system.run("solver", "eval_grad", path=PATH.GRAD,
                   export_traces=PAR.SAVETRACES)
        self.stopwatch("time")

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

    def write_gradient(self):
        """
        Overwrites seisflows.workflow.inversion.write_gradient()

        Removes the need to set path and suffix, hardcodes them into function

        Writes gradient in format expected by nonlinear optimization library
        """
        print("POSTPROCESSING")
        src = os.path.join(PATH.GRAD, "gradient")

        self.stopwatch("set")
        postprocess.write_gradient(PATH.GRAD)
        parts = solver.load(src, suffix="_kernel")
        optimize.save("g_new", solver.merge(parts))
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




