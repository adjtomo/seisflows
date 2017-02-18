
import os
import sys
import numpy as np

from os.path import abspath, basename, join
from seisflows.tools import unix
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class serial(custom_import('system', 'base')):
    """ An interface through which to submit workflows, run tasks in serial or 
      parallel, and perform other system functions.

      By hiding environment details behind a python interface layer, these 
      classes provide a consistent command set across different computing
      environments.

      For important additional information, please see 
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration
    """

    def check(self):
        """ Checks parameters and paths
        """

        # name of job
        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', basename(abspath('.')))

        # number of tasks
        if 'NTASK' not in PAR:
            setattr(PAR, 'NTASK', 1)

        # number of processers per task
        if 'NPROC' not in PAR:
            setattr(PAR, 'NPROC', 1)

        # level of detail in output messages
        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        # where job was submitted
        if 'WORKDIR' not in PATH:
            setattr(PATH, 'WORKDIR', abspath('.'))

        # where output files are written
        if 'OUTPUT' not in PATH:
            setattr(PATH, 'OUTPUT', PATH.WORKDIR+'/'+'output')

        # where temporary files are written
        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH', PATH.WORKDIR+'/'+'scratch')

        # where system files are written
        if 'SYSTEM' not in PATH:
            setattr(PATH, 'SYSTEM', PATH.SCRATCH+'/'+'system')

        # optional local filesystem scratch path
        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)


    def submit(self, workflow):
        """ Submits job
        """
        # create scratch directories
        unix.rm(PATH.SCRATCH)
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # create output directories
        unix.mkdir(PATH.OUTPUT)

        self.checkpoint()

        # execute workflow
        workflow.main()


    def run(self, classname, funcname, hosts='all', **kwargs):
        """ Runs tasks in serial or parallel on specified hosts
        """
        unix.mkdir(PATH.SYSTEM)

        if hosts == 'all':
            for itask in range(PAR.NTASK):
                self.setnode(itask)
                self.progress(itask)
                func = getattr(__import__('seisflows_'+classname), funcname)
                func(**kwargs)
            print ''

        elif hosts == 'head':
            self.setnode(0)
            func = getattr(__import__('seisflows_'+classname), funcname)
            func(**kwargs)

        else:
            task(**kwargs)


    def getnode(self):
        """ Gets number of running task
        """
        return int(os.environ['SEISFLOWS_TASKID'])

    def setnode(self, itask):
        """ Sets number of running task
        """
        os.environ['SEISFLOWS_TASKID'] = str(itask)

    def mpiexec(self):
        """ Specifies MPI exectuable; used to invoke solver
        """
        if PAR.NPROC > 1:
            return 'mpiexec -np %d ' % PAR.NPROC
        else:
            return ''

    def progress(self, itask=None):
        """ Provides status updates
        """
        if PAR.VERBOSE and PAR.NTASK > 1:
            print ' task ' + '%02d'%(itask + 1) + ' of ' + '%02d'%PAR.NTASK
