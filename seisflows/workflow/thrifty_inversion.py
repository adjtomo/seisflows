#!/usr/bin/env python3
"""
Thrifty: using resources carefully and not wastefully

A thrifty inversion skips the costly intialization step (i.e., forward
simulations and misfit quantification) if the final forward simulations from
the previous iteration's line search can be used in the current one.
"""
import sys
import logging

from seisflows.tools import unix, msg
from seisflows.config import custom_import


PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]
optimize = sys.modules["seisflows_optimize"]


class ThriftyInversion(custom_import("workflow", "inversion")):
    """
    Thrifty inversion which attempts to save resources by re-using previous
    line search results for the current iteration.
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        :type thrifty: bool
        :param thrifty: the current status of the inversion.
            if False: assumed to be first iteration, a restart, or some other
            condition has been met which means inversion is defaulting to normal
            behavior
            if True: A well-scaled inversion can skip the function evaluation
            of the next iteration by using the line search results of the
            previous iteration
        """
        super().__init__()
        self.thrifty = False

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=False)
        if validate:
            self.required.validate()

        assert PAR.LINESEARCH == "Backtrack", \
            "Thrifty inversion requires backtracking line search"

    def initialize(self):
        """
        If line search can be carried over, skip initialization step
        Or if manually starting a new run, start with normal inversion init
        """
        if not self.thrifty or optimize.iter == PAR.BEGIN:
            super().initialize()
        else:
            self.logger.info(msg.mjr("INITIALIZING THRIFTY INVERSION"))

    def clean(self):
        """
        Determine if forward simulation from line search can be carried over.
        We assume clean() is the final flow() argument so that we can update
        the thrifty status here.
        """
        self.update_status()

        if self.thrifty:
            self.logger.info(msg.mnr("THRIFTY CLEANING  WORKDIR FOR NEXT "
                                     "ITERATION"))
            unix.rm(PATH.GRAD)
            unix.mv(PATH.FUNC, PATH.GRAD)
            unix.mkdir(PATH.FUNC)
        else:
            super().clean()

    def update_status(self):
        """
        Determine if line search forward simulation can be carried over based
        on a number of criteria
        """
        self.logger.info("updating thrifty inversion status")
        if optimize.iter == PAR.BEGIN:
            self.logger.info("1st iteration, defaulting to inversion workflow")
            thrifty = False
        elif optimize.restarted:
            self.logger.info("optimization has been restarted, defaulting to "
                             "inversion workflow")
            thrifty = False
        elif optimize.iter == PAR.END:
            self.logger.info("final iteration, defaulting to inversion workflow")
            thrifty = False
        else:
            self.logger.info("continuing with thrifty inversion workflow")
            thrifty = True

        self.thrifty = thrifty




