#!/usr/bin/env python3
"""
Thrifty: using resources carefully and not wastefully

A thrifty inversion skips the costly intialization step (i.e., forward
simulations and misfit quantification) if the final forward simulations from
the previous iteration's line search can be used in the current one. Otherwise
it performs the same as the Inversion workflow
"""
import sys
import logging

from seisflows.workflow.inversion import Inversion
from seisflows.tools import unix, msg


class ThriftyInversion(Inversion):
    """
    Thrifty inversion which attempts to save resources by re-using previous
    line search results for the current iteration.
    """
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
        super().check(validate=validate)

        assert self.par.LINESEARCH.upper() == "BACKTRACK", (
            "Thrifty inversion requires PAR.LINESEARCH == 'backtrack'"
        )

    def evaluate_initial_misfit(self):
        """
        If line search can be carried over, skip initialization step
        Or if manually starting a new run, start with normal inversion init
        """
        optimize = self.module("optimize")

        if not self.thrifty or optimize.iter == self.par.BEGIN:
            super().evaluate_initial_misfit()
        else:
            self.logger.info(msg.mjr("SKIPPING INITIAL MISFIT EVALUATION"))

    def clean(self):
        """
        Determine if forward simulation from line search can be carried over.
        We assume clean() is the final flow() argument so that we can update
        the thrifty status here.
        """
        self._update_status()

        if self.thrifty:
            self.logger.info(
                msg.mnr("THRIFTY CLEANING  WORKDIR FOR NEXT ITERATION")
            )
            unix.rm(self.path.GRAD)
            # Last line search evaluation becomes the new gradient evaluation
            unix.mv(self.path.FUNC, self.path.GRAD)
            unix.mkdir(self.path.FUNC)
        else:
            super().clean()

    def _update_status(self):
        """
        Determine if line search forward simulation can be carried over based
        on a variety of criteria relating to location in the inversion.
        """
        optimize = self.module("optimize")
        
        self.logger.info("updating thrifty inversion status")
        if optimize.iter == self.par.BEGIN:
            self.logger.info("1st iteration, defaulting to inversion workflow")
            thrifty = False
        elif optimize.restarted:
            self.logger.info("optimization has been restarted, defaulting to "
                             "inversion workflow")
            thrifty = False
        elif optimize.iter == self.par.END:
            self.logger.info("final iteration, defaulting to inversion workflow")
            thrifty = False
        else:
            self.logger.info("continuing with thrifty inversion workflow")
            thrifty = True

        self.thrifty = thrifty

