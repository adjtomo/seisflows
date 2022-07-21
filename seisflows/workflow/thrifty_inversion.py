#!/usr/bin/env python3
"""
Thrifty: using resources carefully and not wastefully

A thrifty inversion skips the costly intialization step (i.e., forward
simulations and misfit quantification) if the final forward simulations from
the previous iteration's line search can be used in the current one. Otherwise
it performs the same as the Inversion workflow
"""
from seisflows import logger
from seisflows.workflow.inversion import Inversion
from seisflows.tools import unix, msg


class ThriftyInversion(Inversion):
    """
    [workflow.thrifty_inversion] an inversion that attempts to save resources
    by re-using previous line search results for the current iteration.

    :type line_search_method: str
    :param line_search_method: chosen line_search algorithm. Currently available
        are 'bracket' and 'backtrack'. See seisflows.plugins.line_search
        for all available options
    """
    __doc__ = Inversion.__doc__ + __doc__

    def __init__(self, line_search_method):
        """Thrifty does not require input parameters

        """
        super().__init__()

        self._line_search_method = line_search_method
        self._thrifty_status = False

    def check(self):
        """
        Checks that we have the correct line search
        """
        super().check()

        assert(self._line_search_method.title() == "Backtrack"), (
            "Thrifty inversion requires `line_search_method` == 'backtrack'"
        )

    def evaluate_initial_misfit(self):
        """
        If line search can be carried over, skip initialization step
        Or if manually starting a new run, start with normal inversion init
        """
        if not self._thrifty_status or (self.optimize.iteration == self.start):
            super().evaluate_initial_misfit()
        else:
            logger.info(msg.mnr("THRIFTY INVERSION, SKIPPING INITIAL MISFIT "
                                "EVALUATION"))

    def clean_scratch_directory(self):
        """
        Determine if forward simulation from line search can be carried over.
        We assume clean() is the final flow() argument so that we can update
        the thrifty status here.
        """
        self._thrifty_status = self._update_status()

        if self._thrifty_status:
            logger.info(
                msg.mnr("THRIFTY CLEANING  WORKDIR FOR NEXT ITERATION")
            )
            unix.rm(self.path.eval_grad)
            # Last line search evaluation becomes the new gradient evaluation
            unix.mv(self.path.eval_func, self.path.eval_grad)
            unix.mkdir(self.path.eval_func)
        else:
            super().clean_scratch_directory()

    def _update_status(self):
        """
        Determine if line search forward simulation can be carried over based
        on a variety of criteria relating to location in the inversion.
        """
        logger.info("updating thrifty inversion status")
        if self.optimize.iter == self.start:
            logger.info("1st iteration, defaulting to inversion workflow")
            thrifty = False
        elif self.optimize.restarted:
            logger.info("optimization has been restarted, defaulting to "
                        "inversion workflow")
            thrifty = False
        elif self.optimize.iter == self.end:
            logger.info("final iteration, defaulting to inversion workflow")
            thrifty = False
        else:
            logger.info("continuing with thrifty inversion workflow")
            thrifty = True

        return thrifty

