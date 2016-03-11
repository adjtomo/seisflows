
from os.path import join
import sys
import numpy as np

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import divides, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import optimize
import preprocess
import postprocess


class thrifty_inversion(custom_import('workflow', 'inversion')):
    """ Thrifty inversion subclass

      Provides savings over conventional inversion by avoiding redundant forward
      simulations associated with sufficient decrease and curvature tests in a 
      safeguarded backtracking line search.

      The results of 'inversion' and 'thrifty_inversion' should be exactly the
      same.  Users who prefer a simpler but less efficient workflow can choose
      choose 'inversion'. Users who prefer a more efficient but more complicated
      workflow can choose 'thrifty_inversion.'
    """

    def solver_status(self, maxiter=1):
        """ Keeps track of whether a forward simulation would be redundant
        """
        if optimize.iter <= maxiter:
            # forward simulation not redundant because solver files do not exist
            # prior to first iteration
            return False

        elif optimize.iter == PAR.BEGIN:
            # forward simulation not redundant because solver files need to be
            # reinstated after possible multiscale transition
            return False

        elif PATH.LOCAL:
            # forward simulation not redundant because solver files need to be
            # reinstated on local filesystems
            return False

        elif PAR.LINESEARCH != 'Backtrack':
            # thrifty inversion only implemented for backtracking line search,
            # not bracketing line search
            return False

        elif optimize.restarted:
            # forward simulation not redundant following optimization algorithm
            # restart
            return False

        else:
            # if none of the above conditions are triggered, then forward 
            # simulation is redundant, can be skipped
            return True

    def setup(self):
        """ Lays groundwork for inversion
        """
        # clean scratch directories
        if PAR.BEGIN == 1:
            unix.rm(PATH.SCRATCH)
            unix.mkdir(PATH.SCRATCH)

            preprocess.setup()
            postprocess.setup()
            optimize.setup()

        isready = self.solver_status()
        if not isready:
            if PATH.DATA:
                print 'Copying data'
            else:
                print 'Generating data'

            system.run('solver', 'setup',
                       hosts='all')


    def initialize(self):
        # are prerequisites for gradient evaluation in place?
        isready = self.solver_status(maxiter=2)

        # if not, then prepare for gradient evaluation
        if not isready:
            super(thrifty_inversion, self).initialize()


    def iterate_search(self):
        super(thrifty_inversion, self).iterate_search()

        isdone = optimize.isdone
        isready = self.solver_status()

        # to avoid redundant forward simulation, save solver files associated
        # with 'best' trial model
        if isready and isdone:
            unix.rm(PATH.SOLVER+'_best')
            unix.mv(PATH.SOLVER, PATH.SOLVER+'_best')


    def clean(self):
        isready = self.solver_status()
        if isready:
            unix.rm(PATH.GRAD)
            unix.mv(PATH.FUNC, PATH.GRAD)
            unix.mkdir(PATH.FUNC)
            unix.rm(PATH.SOLVER)
            unix.mv(PATH.SOLVER+'_best', PATH.SOLVER)
        else:
            super(thrifty_inversion, self).clean()


