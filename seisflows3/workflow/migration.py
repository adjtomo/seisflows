#!/usr/bin/env python
"""
This is the base class seisflows.workflow.migration

This is a main Seisflows class, it controls the main workflow.
"""
import os
import sys
import logging

from seisflows3.tools import unix, msg
from seisflows3.tools.wrappers import exists
from seisflows3.config import custom_import, SeisFlowsPathsParameters


PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]

system = sys.modules["seisflows_system"]
solver = sys.modules["seisflows_solver"]
preprocess = sys.modules["seisflows_preprocess"]
postprocess = sys.modules["seisflows_postprocess"]


class Migration(custom_import("workflow", "base")):
    """
    Migration base class.

    Performs the workflow of an inversion up to the postprocessing. In the
    terminology of seismic exploration, implements a 'reverse time migration'.
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        """
        super().__init__()

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters(super().required)

        return sf

    def main(self, return_flow=False):
        """
        Migrates seismic data
        """
        flow = [self.setup,
                self.generate_synthetics,
                self.backproject,
                self.process_kernels
                ]
        if return_flow:
            return flow
        else:
            # FLOW is a constant, when we run, we flow
            flow = FLOW

        # Run each argument in flow
        for func in flow:
            func()

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVEKERNELS:
            self.save_kernels()
        else:
            self.save_kernels_sum()

    def setup(self):
        """
        Sets up the SeisFlows3 modules for the Migration
        """
        # Set up all the requisite modules from the master job
        self.logger.info(msg.mnr("PERFORMING MODULE SETUP"))
        preprocess.setup()
        postprocess.setup()
        system.run("solver", "setup")

    def generate_synthetics(self):
        """
        Performs forward simulation, and evaluates the objective function
        """
        self.logger.info(msg.sub("PREPARING VELOCITY MODEL"))
        src = os.path.join(PATH.OUTPUT, "model_init")
        dst = os.path.join(PATH.SCRATCH, "model")

        assert os.path.exists(src)
        unix.cp(src, dst)

        self.logger.info(msg.sub("EVALUATE OBJECTIVE FUNCTION"))
        system.run("solver", "eval_func", path=PATH.SCRATCH,
                   write_residuals=True)

    def backproject(self):
        """
        Backproject or create kernels by running adjoint simulations
        """
        self.logger.info(msg.sub("BACKPROJECT / EVALUATE GRADIENT"))
        system.run("solver", "eval_grad", path=PATH.SCRATCH,
                   export_traces=PAR.SAVETRACES)

    def process_kernels(self):
        """
        Backproject to create kernels from synthetics
        """
        system.run_single("postprocess", "process_kernels",
                          path=os.path.join(PATH.SCRATCH, "kernels"),
                          parameters=solver.parameters)

        try:
            # TODO Figure out a better method for running this try except
            system.run_single("postprocess", "process_kernels",
                              path=os.path.join(PATH.SCRATCH, "kernels"),
                              parameters=["rhop"])
        except:
            pass

    def save_kernels_sum(self):
        src = PATH.SCRATCH +"/"+ "kernels/sum"
        dst = PATH.OUTPUT +"/"+ "kernels"
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_kernels(self):
        src = PATH.SCRATCH +"/"+ "kernels"
        dst = PATH.OUTPUT
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_traces(self):
        src = PATH.SCRATCH +"/"+ "traces"
        dst = PATH.OUTPUT
        unix.cp(src, dst)

