#!/usr/bin/env python
"""
This is the Base class for seisflows.workflow.
It contains mandatory functions that must be called by subclasses
"""
import sys
import logging

from seisflows3.tools import msg
from seisflows3.tools.wrappers import exists
from seisflows3.config import save, SeisFlowsPathsParameters


PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]


class Base:
    """
    Workflow abstract base class
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by the user.
        Attributes are initialized as NoneTypes for clarity and docstrings.
        """
        pass

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        sf.par("CASE", required=True, par_type=str,
               docstr="Type of inversion, available: "
                      "['data': real data inversion, "
                      "'synthetic': synthetic-synthetic inversion]")

        sf.par("RESUME_FROM", required=False, par_type=str,
               docstr="Name of task to resume inversion from")

        sf.par("STOP_AFTER", required=False, par_type=str,
               docstr="Name of task to stop inversion after finishing")

        sf.par("SAVEMODEL", required=False, default=True, par_type=bool,
               docstr="Save final model files after each iteration")

        sf.par("SAVEGRADIENT", required=False, default=True, par_type=bool,
               docstr="Save gradient files after each iteration")

        sf.par("SAVEKERNELS", required=False, default=False, par_type=bool,
               docstr="Save event kernel files after each iteration")

        sf.par("SAVETRACES", required=False, default=False, par_type=bool,
               docstr="Save waveform traces after each iteration")

        sf.par("SAVERESIDUALS", required=False, default=False, par_type=bool,
               docstr="Save waveform residuals after each iteration")

        sf.par("SAVEAS", required=False, default="binary", par_type=str,
               docstr="Format to save models, gradients, kernels. "
                      "Available: "
                      "['binary': save files in native SPECFEM .bin format, "
                      "'vector': save files as NumPy .npy files, "
                      "'both': save as both binary and vectors]")

        sf.path("MODEL_INIT", required=True,
                docstr="location of the initial model to be used for workflow")

        sf.path("MODEL_TRUE", required=False,
                docstr="Target model to be used for PAR.CASE == 'synthetic'")

        sf.path("DATA", required=False, default=None,
                docstr="path to data available to workflow")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths. Must be implemented by sub-class
        """
        msg.check(type(self))
        if validate:
            self.required.validate()

        if PAR.CASE.upper() == "SYNTHETIC":
            assert exists(PATH.MODEL_TRUE), \
                "CASE == SYNTHETIC requires PATH.MODEL_TRUE"

        if not exists(PATH.DATA):
            assert "MODEL_TRUE" in PATH, f"DATA or MODEL_TRUE must exist"

    def main(self, flow, return_flow=False):
        """
        Execution of a workflow is equal to stepping through workflow.main()
        """
        pass

    def resume_from(self, flow):
        """
        Allow the main() function to resume a workflow from a given flow
        argument. Allows User to attempt restarting failed or stopped workflow
        """
        try:
            resume_idx = [_.__name__ for _ in flow].index(PAR.RESUME_FROM)
            resume_fx = flow[resume_idx].__name__
        except ValueError:
            self.logger.info(f"{PAR.RESUME_FROM} does not correspond to any "
                             f"workflow functions. Exiting...")
            sys.exit(-1)

        self.logger.info(
            msg.mjr(f"RESUME WORKFLOW FROM FUNCTION: '{resume_fx}'")
        )
        # Curtail the flow argument during resume_from
        return flow[resume_idx:]

    @staticmethod
    def checkpoint():
        """
        Writes information to disk so workflow can be resumed following a break
        """
        save()


