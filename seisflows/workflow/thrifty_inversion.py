
import sys

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.config import ParameterError, custom_import
from seisflows.workflow.base import base

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

optimize = sys.modules['seisflows_optimize']


class thrifty_inversion(custom_import('workflow', 'inversion')):
    """ Thrifty inversion subclass

      Provides savings over conventional inversion by carrying over forward
      simulations from line search

      The results of 'inversion' and 'thrifty_inversion' should be exactly the
      same
    """

    status=0

    def initialize(self):
        if self.status==0:
            super(thrifty_inversion, self).initialize()


    def clean(self):
        # can forward simulations from line search be carried over?
        self.update_status()

        if self.status==1:
            unix.rm(PATH.GRAD)
            unix.mv(PATH.FUNC, PATH.GRAD)
            unix.mkdir(PATH.FUNC)
        else:
            super(thrifty_inversion, self).clean()


    def update_status(self):
        if PAR.LINESEARCH != 'Backtrack':
            # only works for backtracking line search
            self.status=0

        elif optimize.iter==PAR.BEGIN or \
             optimize.restarted:
            # even if backtracking line search is chosen, may not work on
            # first iteration or following a restart
            self.status=0

        elif optimize.iter==PAR.END:
            # may not work after resuming saved workflow
            self.status=0

        elif PATH.LOCAL:
            # may not work if using local filesystems
            self.status=0

        else:
            self.status=1



