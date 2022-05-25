#!/usr/bin/env python3
"""
This is the Base class for seisflows.workflow.
It contains mandatory functions that must be called by subclasses
"""
import sys
import logging

from seisflows.tools import msg
from seisflows.tools.wrappers import exists
from seisflows.config import save, SeisFlowsPathsParameters


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
        if validate:
            self.required.validate()

        if PAR.CASE.upper() == "SYNTHETIC":
            assert exists(PATH.MODEL_TRUE), \
                "CASE == SYNTHETIC requires PATH.MODEL_TRUE"

        if not exists(PATH.DATA):
            assert "MODEL_TRUE" in PATH, f"DATA or MODEL_TRUE must exist"

    def main(self, return_flow=False):
        """
        Execution of a workflow is equal to stepping through workflow.main()

        An example main() script is provided below which details the requisite
        parts. This function will NOT execute as it is written in pseudocode.

        :type return_flow: bool
        :param return_flow: for CLI tool, simply returns the flow function
            rather than running the workflow. Used for print statements etc.
        """
        self.logger.warning("The current definition of workflow.main() will "
                            "NOT execute, it must be overwritten by a "
                            "subclass")

        # The FLOW function defines a list of functions to execute IN ORDER
        flow = (self.func1,
                self.func2,
                # ...
                self.funcN
                )

        # REQUIRED: CLI command `seisflows print flow` needs this for output
        if return_flow:
            return flow

        # Allow User to start the workflow mid-FLOW, in the event that a
        # previous workflow errored, or if the User had previously stopped
        # a workflow to look at results and they want to pick up where
        # they left off
        start, stop = self.check_stop_resume_cond(flow)

        self.logger.info(msg.mjr("BEGINNING EXAMPLE WORKFLOW"))

        # Iterate through the `FLOW` to step through workflow.main()
        for func in flow[start:stop]:
            func()

            # If an multi-iteration workflow is run, `FLOW` will be executed
            # repeatedly, so reset start and stop for subsequent iterations
            start, stop = 0, -1

        self.logger.info(msg.mjr("FINISHED EXAMPLE WORKFLOW"))

    def check_stop_resume_cond(self, flow):
        """
        Chek the stop after and resume from conditions

        Allow the main() function to resume a workflow from a given flow
        argument, or stop the workflow after a given argument. In the event
        that a previous workflow errored, or if the User had previously
        stopped a workflow to look at results and they want to pick up where
        they left off.

        Late check: Exits the workflow if RESUME_FROM or STOP_AFTER arguments
        do not match any of the given flow arguments.

        :type flow: tuple of functions
        :param flow: an ordered list of functions that will be
        :rtype: tuple of int
        :return: (start, stop) indices of the `flow` input dictating where the
            list should be begun and ended. If RESUME_FROM and STOP_AFTER
            conditions are NOT given by the user, start and stop will be 0 and
            -1 respectively, meaning we should execute the ENTIRE list
        """
        fxnames = [func.__name__ for func in flow]

        # Default values which dictate that flow will execute in its entirety
        start_idx = None
        stop_idx = None

        # Overwrite start_idx if RESUME_FROM given, exit condition if no match
        if PAR.RESUME_FROM:
            try:
                start_idx = fxnames.index(PAR.RESUME_FROM)
                fx_name = flow[start_idx].__name__
                self.logger.info(
                    msg.mnr(f"WORKFLOW WILL RESUME FROM FUNC: '{fx_name}'")
                )
            except ValueError:
                self.logger.info(
                    msg.cli(f"{PAR.RESUME_FROM} does not correspond to any FLOW "
                            f"functions. Please check that PAR.RESUME_FROM "
                            f"matches one of the functions listed out in "
                            f"`seisflows print flow`.", header="error",
                            border="=")
                )
                sys.exit(-1)

        # Overwrite stop_idx if STOP_AFTER provided, exit condition if no match
        if PAR.STOP_AFTER:
            try:
                stop_idx = fxnames.index(PAR.STOP_AFTER)
                fx_name = flow[stop_idx].__name__
                stop_idx += 1  # increment to stop AFTER, due to python indexing
                self.logger.info(
                    msg.mnr(f"WORKFLOW WILL STOP AFTER FUNC: '{fx_name}'")
                )
            except ValueError:
                self.logger.info(
                    msg.cli(f"{PAR.STOP_AFTER} does not correspond to any FLOW "
                            f"functions. Please check that PAR.STOP_AFTER "
                            f"matches one of the functions listed out in "
                            f"`seisflows print flow`.", header="error",
                            border="=")
                )
                sys.exit(-1)

        # Make sure stop after doesn't come before resume_from, otherwise none
        # of the flow will execute
        if PAR.STOP_AFTER and PAR.RESUME_FROM:
            if stop_idx <= start_idx:
                self.logger.info(
                    msg.cli(f"PAR.STOP_AFTER=='{PAR.STOP_AFTER}' is called "
                            f"before PAR.RESUME_FROM=='{PAR.RESUME_FROM}' in "
                            f"the FLOW functions. Please adjust accordingly "
                            f"and rerun.", header="error", border="=")
                )
                sys.exit(-1)

        return start_idx, stop_idx

    @staticmethod
    def checkpoint():
        """
        Writes information to disk so workflow can be resumed following a break
        """
        save()


