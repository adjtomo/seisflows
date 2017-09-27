
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

        # how to invoke executables
        if 'MPIEXEC' not in PAR:
            setattr(PAR, 'MPIEXEC', '')

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
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # create output directories
        unix.mkdir(PATH.OUTPUT)

        workflow.checkpoint()

        # execute workflow
        workflow.main()


    def run(self, classname, method, hosts='all', **kwargs):
        """ Executes task multiple times in serial
        """
        unix.mkdir(PATH.SYSTEM)

        for taskid in range(PAR.NTASK):
            os.environ['SEISFLOWS_TASKID'] = str(taskid)
            if PAR.VERBOSE > 0:
                self.progress(taskid)
            func = getattr(__import__('seisflows_'+classname), method)
            func(**kwargs)
        print ''


    def run_single(self, classname, method, *args, **kwargs):
        """ Runs task a single time
        """
        os.environ['SEISFLOWS_TASKID'] = str(0)
        func = getattr(__import__('seisflows_'+classname), method)
        func(**kwargs)


    def taskid(self):
        """ Provides a unique identifier for each running task
        """
        return int(os.environ['SEISFLOWS_TASKID'])


    def mpiexec(self):
        """ Specifies MPI executable used to invoke solver
        """
        return PAR.MPIEXEC


    def progress(self, taskid):
        """ Provides status update
        """
        if PAR.NTASK > 1:
            print ' task ' + '%02d of %02d' % (taskid+1, PAR.NTASK)
