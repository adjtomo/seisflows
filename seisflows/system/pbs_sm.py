
import os
import math
import sys
import time

from os.path import abspath, basename, join
from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.tools import call, findpath, saveobj
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class pbs_lg(custom_import('system', 'base')):
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
        print msg.Warning_pbs_sm

        # name of job
        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', basename(abspath('.')))

        # time allocated for workflow in minutes
        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        # number of tasks
        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        # number of cores per task
        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        # number of cores per node
        if 'NODESIZE' not in PAR:
            raise ParameterError(PAR, 'NODESIZE')

        # how to invoke executables
        if 'MPIEXEC' not in PAR:
            setattr(PAR, 'MPIEXEC', '')

        # optional additional PBS arguments
        if 'PBSARGS' not in PAR:
            setattr(PAR, 'PBSARGS', '')

        # optional environment variable list VAR1=val1,VAR2=val2,...
        if 'ENVIRONS' not in PAR:
            setattr(PAR, 'ENVIRONS', '')

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

        # optional local scratch path
        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)


    def submit(self, workflow):
        """ Submits workflow
        """
        # create scratch directories
        unix.mkdir(PATH.SCRATCH)
        unix.mkdir(PATH.SYSTEM)

        # create output directories
        unix.mkdir(PATH.OUTPUT)

        workflow.checkpoint()

        hours = PAR.WALLTIME/60
        minutes = PAR.WALLTIME%60
        walltime = 'walltime=%02d:%02d:00 ' % (hours, minutes)

        ncpus = PAR.NODESIZE
        mpiprocs = PAR.NODESIZE

        # prepare qsub arguments
        call( 'qsub '
                + '%s ' % PAR.PBSARGS
                + '-l select=1:ncpus=%d:mpiprocs=%d ' % (ncpus, mpiprocs)
                + '-l %s ' % walltime
                + '-N %s ' % PAR.TITLE
                + '-j %s '%'oe'
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.log')
                + '-V '
                + ' -- ' + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


    def run(self, classname, method, hosts='all', **kwargs):
        """ Runs embarrassingly parallel tasks

          Executes the following multiple times:
              classname.method(*args, **kwargs)

          system.taskid serves to provide each running task a unique identifier
        """
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        if hosts == 'all':
            # run all tasks
            call(findpath('seisflows.system')  +'/'+'wrappers/dsh '
                    + ','.join(self.hostlist()) + ' '
                    + findpath('seisflows.system')  +'/'+'wrappers/run '
                    + PATH.OUTPUT + ' '
                    + classname + ' '
                    + method + ' '
                    + 'PYTHONPATH='+findpath('seisflows'),+','
                    + PAR.ENVIRONS)

        elif hosts == 'head':
            # run a single task
            call('ssh ' + self.hostlist()[0] + ' '
                    + '"'
                    + 'export SEISFLOWS_TASK_ID=0; '
                    + join(findpath('seisflows.system'), 'wrappers/run ')
                    + PATH.OUTPUT + ' '
                    + classname + ' '
                    + method + ' '
                    + 'PYTHONPATH='+findpath('seisflows'),+','
                    + PAR.ENVIRONS
                    +'"')

        else:
            raise KeyError('Bad keyword argument: system.run: hosts')



    def mpiexec(self):
        """ Specifies MPI executable used to invoke solver
        """
        return PAR.MPIEXEC


    def taskid(self):
        """ Provides a unique identifier for each running task
        """
        try:
            return os.getenv('PBS_NODENUM')
        except:
            raise Exception("PBS_NODENUM environment variable not defined.")


    def hostlist(self):
        """ Generates list of allocated cores
        """
        with open(os.environ['PBS_NODEFILE'], 'r') as f:
            return [line.strip() for line in f.readlines()]


    def save_kwargs(self, classname, method, kwargs):
        kwargspath = join(PATH.OUTPUT, 'kwargs')
        kwargsfile = join(kwargspath, classname+'_'+method+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)


