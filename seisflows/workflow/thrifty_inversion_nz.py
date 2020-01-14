
import sys

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.config import ParameterError, custom_import
from seisflows.workflow.base import base

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

optimize = sys.modules['seisflows_optimize']


class thrifty_inversion_nz(custom_import('workflow', 'inversion_nz')):
    """ Thrifty inversion subclass

      Provides savings over conventional inversion by carrying over forward
      simulations from line search

      The results of 'inversion' and 'thrifty_inversion' should be exactly the
      same
    """
    # Instanstiate the status attribute
    status = 0

    def initialize(self):
        """
        If line search can be carried over, skip initialization step
        Or if manually starting a new run, start with normal inversion init
        """
        if (self.status == 0) or (optimize.iter == PAR.BEGIN):
            super(thrifty_inversion_nz, self).initialize()
        else:
            print 'THRIFTY INITIALIZE'


    def clean(self):
        """
        Determine if forward simulation from line search can be carried over
        """
        self.update_status()
        if self.status == 1:
            print 'THRIFTY CLEAN'
            unix.rm(PATH.GRAD)
            unix.mv(PATH.FUNC, PATH.GRAD)
            unix.mkdir(PATH.FUNC)
        else:
            super(thrifty_inversion_nz, self).clean()


    def update_status(self):
        """
        Determine if line search forward simulation can be carried over
        """
        print 'THRIFTY STATUS'
        # only works for backtracking line search
        if PAR.LINESEARCH != 'Backtrack':
            print '\t Line search not "Backtrack", cannot run thrifty'
            self.status = 0
        # may not work on first iteration
        elif optimize.iter == PAR.BEGIN:
            print '\t First iteration of workflow, defaulting to inversion'
            self.status=0
        # may not work following restart
        elif optimize.restarted:
            print '\t Optimization has been restarted, defaulting to inversion'
            self.status = 0
        # may not work after resuming saved workflow
        elif optimize.iter == PAR.END:
            print '\t End of workflow, defaulting to inversion'
            self.status = 0
        # may not work if using local filesystems
        elif PATH.LOCAL:
            print '\t Local filesystem, cannot run thrifty'
            self.status = 0
        # otherwise, continue with thrifty inversion
        else:
            print '\t Continuing with thrifty inversion'
            self.status = 1



