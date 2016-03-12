
import os
import math
import sys
import time
from os.path import abspath, basename, join

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.code import call, findpath, saveobj
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class lsf_sm(custom_import('system', 'mpi')):
    """ An interface through which to submit workflows, run tasks in serial or 
      parallel, and perform other system functions.

      By hiding environment details behind a python interface layer, these 
      classes provide a consistent command set across different computing
      environments.

      Intermediate files are written to a global scratch path PATH.SCRATCH,
      which must be accessible to all compute nodes.

      Optionally, users can provide a local scratch path PATH.LOCAL if each
      compute node has its own local filesystem.

      For important additional information, please see 
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(lsf_sm, self).check()

        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        if 'LSFARGS' not in PAR:
            setattr(PAR, 'LSFARGS', '')


    def submit(self, workflow):
        """ Submits workflow
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)
        unix.mkdir(PATH.SUBMIT+'/'+'output.lsf')

        self.checkpoint()

        # prepare bsub arguments
        call('bsub '
                + '%s ' % PAR.LSFARGS
                + '-J %s ' % PAR.TITLE
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.log')
                + '-n %d ' % PAR.NTASK
                + '-e %s ' % (PATH.SUBMIT+'/'+'error.log')
                + '-W %d:00 ' % PAR.WALLTIME
                + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


