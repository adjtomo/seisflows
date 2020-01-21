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

from pyatoa.plugins.seisflows.pyaflowa import Pyaflowa

# Seisflows configuration
PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']
optimize = sys.modules['seisflows_optimize']
preprocess = sys.modules['seisflows_preprocess']
postprocess = sys.modules['seisflows_postprocess']


class InversionPyatoa(custom_import('workflow', 'Inversion')):
    """
    Waveform inversion subclass, additional support for Pyatoa integration
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

            print("\tInitializing Pyatoa")
            pyaflowa = Pyaflowa(seisflows_paths=PATH,
                                seisflows_parameters=PAR)
            pyaflowa.setup()
            
            print("\tPreparing initial model")
            self._watch("set")
            system.run("solver", "setup")
            self._watch("time")

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

        print("\tRunning forward simulation")
        self._watch("set")
        system.run("solver", "eval_fwd", path=path_)
        self._watch("time")

        print("\tQuantifying misfit")
        self._watch("set")
        self.pyaflowa.set(iter=optimize.iter, suffix=suffix_)
        system.run_ancil("solver", "eval_func", pyaflowa=self.pyaflowa)
        self._watch("time")

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

        print("\tRunning forward simulation")
        self._watch("set")
        system.run("solver", "eval_fwd", path=path_)
        self._watch("time")

        print("\tQuantifying misfit")
        self._watch("set")
        self.pyaflowa.set(iter=optimize.iter,
                          step=optimize.line_search.step_count + 1,
                          suffix=suffix_)
        system.run_ancil("solver", "eval_func", pyaflowa=self.pyaflowa)
        self._watch("time")

        self.write_misfit(suffix=suffix_)

    def finalize(self):
        """
        Overwrites seisflows.workflow.inversion.finalize()

        Adds finalization of Pyaflowa and saves results from current iteration
        """
        super(InversionPyatoa, self).finalize()

        # Finalize Pyatoa run
        self.pyaflowa.set(iter=optimize.iter,
                          model=f"m{int(optimize.iter)-1:0>2}",
                          step=f"s{optimize.line_search.step_count:0>2}")
        self.pyaflowa.finalize()

    def write_misfit(self, suffix=''):
        """
        Overwrites seisflows.workflow.inversion.write_misfit()

        Writes misfit in format expected by nonlinear optimization library
        Waits for all instances of Pyatoa to finish writing their misfit

        !!! This should be rewritten, it's a bit clunky !!!

        :type suffix: str
        :param suffix: suffix to add to the misfit
        """
        src = os.path.join(PATH.PYATOA_IO, "data", "misfits", "*")
        dst = "f_{suffix}"

        # Used to wait for processes to finish
        misfit = 0
        waited_s = 0
        wait_interval_s = 5

        while True:
            misfits = glob.glob(src)
            if len(misfits) == PAR.NSRC:
                for _misfit in misfits:
                    misfit += np.loadtxt(_misfit)
                    os.remove(_misfit)
                # Following Tape et al 2007, total misfit=misfit/number_srcs
                optimize.savetxt(dst, misfit/PAR.NSRC)
                return
            else:
                waited_s += wait_interval_s
                # If this function has waited more than x minutes, exit because
                # something has gone wrong with misfit quantification
                if waited_s >= (60 * 5):
                    print("Error, misfit writer stuck waiting for files")
                    sys.exit(-1)

                time.sleep(wait_interval_s)   





