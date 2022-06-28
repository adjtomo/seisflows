#!/usr/bin/env python3
"""
This is the base class seisflows.workflow.migration

This is a main Seisflows class, it controls the main workflow.
"""
import os

from seisflows.workflow.forward import Forward
from seisflows.tools import unix, msg


class Migration(Forward):
    """
    Migration base class.

    Performs the workflow of an inversion up to the postprocessing. In the
    terminology of seismic exploration, implements a 'reverse time migration'.
    """
    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.

        """
        super().__init__()

    def check(self, validate=True):
        """
        Checks parameters and paths. Must be implemented by sub-class
        """
        super().check(validate=validate)

    def setup(self, flow=None, return_flow=False):
        """
        Override the Forward.setup() method to include new flow functions
        AND run setup for a the Postprocess module which will be used to deal
        with the gradient
        """
        if flow is None:
            flow = (self.evaluate_initial_misfit,
                    self.evaluate_gradient)

        super().setup(flow=flow, return_flow=return_flow)

        postprocess = self.module("postprocess")
        postprocess.setup()

    def main(self, flow=None, return_flow=False):
        """Inherits from seisflows.workflow.forward.Forward"""
        self.main(flow=flow, return_flow=return_flow)

    def evaluate_initial_misfit(self):
        """Inherits from seisflows.workflow.forward.Forward"""
        self.evaluate_initial_misfit()

    def evaluate_gradient(self, path=None):
        """
        Performs adjoint simulation to retrieve the gradient of the objective
        """
        system = self.module("system")
        
        self.logger.info(msg.mnr("EVALUATING GRADIENT"))

        self.logger.debug(f"evaluating gradient {self.par.NTASK} times on system...")
        system.run("solver", "eval_grad", path=path or self.path.GRAD,
                   export_traces=self.par.SAVETRACES)

    def process_kernels(self):
        """
        Backproject to create kernels from synthetics
        """
        system = self.module("system")
        solver = self.module("solver")

        system.run("postprocess", "process_kernels", single=True,
                   path=os.path.join(self.path.SCRATCH, "kernels"),
                   parameters=solver.parameters)

        try:
            # TODO Figure out a better method for running this try except
            system.run("postprocess", "process_kernels", single=True,
                       path=os.path.join(self.path.SCRATCH, "kernels"),
                       parameters=["rhop"])
        except:
            pass

    def finalize(self):
        """
        Saves results from current model update iteration
        """
        self.logger.info(msg.mnr("FINALIZING MIGRATION WORKFLOW"))

        if self.par.SAVETRACES:
            self.save_traces()
        if self.par.SAVEKERNELS:
            self.save_kernels()
        else:
            self.save_kernels_sum()

    def save_kernels_sum(self):
        """
        Same summed kernels into the output directory
        """
        src = os.path.join(self.path.SCRATCH, "kernels", "sum")
        dst = os.path.join(self.path.OUTPUT, "kernels")
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_kernels(self):
        """
        Save individual kernels into the output directory
        """
        src = os.path.join(self.path.SCRATCH, "kernels")
        dst = self.path.OUTPUT
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_traces(self):
        """
        Save waveform traces into the output directory
        """
        src = os.path.join(self.path.SCRATCH, "traces")
        dst = self.path.OUTPUT
        unix.cp(src, dst)

