#!/usr/bin/env python3
"""
This is the base class seisflows.workflow.migration

This is a main Seisflows class, it controls the main workflow.
"""
import os
import sys
import logging

from seisflows.tools import unix, msg
from seisflows.tools.wrappers import exists
from seisflows.config import custom_import, SeisFlowsPathsParameters


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
        """s
        Migrates seismic data to generate sensitivity kernels

        :type return_flow: bool
        :param return_flow: for CLI tool, simply returns the flow function
            rather than running the workflow. Used for print statements etc.
        """
        flow = (self.setup,
                self.generate_synthetics,
                self.backproject,
                self.process_kernels,
                self.finalize,
                )
        if return_flow:
            return flow

        # Allow workflow resume from and stop after given flow functions
        start, stop = self.check_stop_resume_cond(flow)

        # Run each argument in flow
        self.logger.info(msg.mjr("STARTING MIGRATION WORKFLOW"))
        for func in flow[start:stop]:
            func()
        self.logger.info(msg.mjr("FINISHED MIGRATION WORKFLOW"))

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
        system.run("postprocess", "process_kernels", single=True,
                   path=os.path.join(PATH.SCRATCH, "kernels"),
                   parameters=solver.parameters)

        try:
            # TODO Figure out a better method for running this try except
            system.run("postprocess", "process_kernels", single=True,
                       path=os.path.join(PATH.SCRATCH, "kernels"),
                       parameters=["rhop"])
        except:
            pass

    def finalize(self):
        """
        Saves results from current model update iteration
        """
        self.logger.info(msg.mnr("FINALIZING MIGRATION WORKFLOW"))

        if PAR.SAVETRACES:
            self.save_traces()
        if PAR.SAVEKERNELS:
            self.save_kernels()
        else:
            self.save_kernels_sum()

    def save_kernels_sum(self):
        """
        Same summed kernels into the output directory
        """
        src = os.path.join(PATH.SCRATCH, "kernels", "sum")
        dst = os.path.join(PATH.OUTPUT, "kernels")
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_kernels(self):
        """
        Save individual kernels into the output directory
        """
        src = os.path.join(PATH.SCRATCH, "kernels")
        dst = PATH.OUTPUT
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_traces(self):
        """
        Save waveform traces into the output directory
        """
        src = os.path.join(PATH.SCRATCH, "traces")
        dst = PATH.OUTPUT
        unix.cp(src, dst)

