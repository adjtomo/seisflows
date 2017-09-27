
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
        print msg.Warning_pbs_lg

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
            setattr(PAR, 'MPIEXEC', 'mpiexec')

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
        unix.mkdir(PATH.WORKDIR+'/'+'output.pbs')

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
                + '-o %s ' % (PATH.WORKDIR+'/'+'output.log')
                + '-V '
                + ' -- ' + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


    def run(self, classname, method, hosts='all', **kwargs):
        """ Runs task multiple times in embarrassingly parallel fasion

          Executes classname.method(*args, **kwargs) NTASK times, each time on
          NPROC cpu cores
        """
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        jobs = self.submit_job_array(classname, method, hosts)
        while True:
            # wait a few seconds before checking again
            time.sleep(5)
            self._timestamp()
            isdone, jobs = self.job_array_status(classname, method, jobs)
            if isdone:
                return


    def mpiexec(self):
        """ Specifies MPI executable used to invoke solver
        """
        return PAR.MPIEXEC


    def taskid(self):
        """ Provides a unique identifier for each running task
        """
        try:
            return os.getenv('PBS_ARRAY_INDEX')
        except:
            raise Exception("PBS_ARRAY_INDEX environment variable not defined.")


    ### private methods

    def submit_job_array(self, classname, method, hosts='all'):
        with open(PATH.SYSTEM+'/'+'job_id', 'w') as f:
            call(self.job_array_cmd(classname, method, hosts),
                stdout=f)

        # retrieve job ids
        with open(PATH.SYSTEM+'/'+'job_id', 'r') as f:
            line = f.readline()
            job = line.split()[-1].strip()
        if hosts == 'all' and PAR.NTASK > 1:
            nn = range(PAR.NTASK)
            job0 = job.strip('[].sdb')
            return [job0+'['+str(ii)+'].sdb' for ii in nn]
        else:
            return [job]


    def job_array_cmd(self, classname, method, hosts):
        nodes = math.ceil(PAR.NTASK/float(PAR.NODESIZE))
        ncpus = PAR.NPROC
        mpiprocs = PAR.NPROC

        hours = PAR.TASKTIME/60
        minutes = PAR.TASKTIME%60
        walltime = 'walltime=%02d:%02d:00 '%(hours, minutes)

        return ('qsub '
                + '%s ' % PAR.PBSARGS
                + '-l select=%d:ncpus=%d:mpiprocs=%d ' (nodes, ncpus, mpiprocs)
                + '-l %s ' % walltime
                + '-J 0-%s ' % (PAR.NTASK-1)
                + '-N %s ' % PAR.TITLE
                + '-o %s ' % (PATH.WORKDIR+'/'+'output.pbs/' + '$PBS_ARRAYID')
                + '-r y '
                + '-j oe '
                + '-V '
                + self.job_array_args(hosts)
                + PATH.OUTPUT + ' '
                + classname + ' '
                + method + ' '
                + 'PYTHONPATH='+findpath('seisflows.system'),+','
                + PAR.ENVIRONS)


    def job_array_args(self, hosts):
        if hosts == 'all':
          args = ('-J 0-%s ' % (PAR.NTASK-1)
                +'-o %s ' % (PATH.WORKDIR+'/'+'output.pbs/' + '$PBS_ARRAYID')
                + ' -- ' + findpath('seisflows.system') +'/'+ 'wrappers/run ')

        elif hosts == 'head':
          args = ('-J 0-0 '
                 +'-o %s ' % (PATH.WORKDIR+'/'+'output.pbs/' + '$PBS_JOBID')
                 + ' -- ' + findpath('seisflows.system') +'/'+ 'wrappers/run ')

        return args


    def job_array_status(self, classname, method, jobs):
        """ Determines completion status of one or more jobs
        """
        states = []
        for job in jobs:
            state = self._query(job)
            if state in ['C']:
                states += [1]
            else:
                states += [0]
            if state in ['F']:
                print msg.TaskError_PBS % (classname, method, job)
                sys.exit(-1)
        isdone = all(states)

        return isdone, jobs


    def _query(self, jobid):
        """ Queries job state from PBS database
        """
        # TODO: replace shell utilities with native Python
        with open(PATH.SYSTEM+'/'+'job_status', 'w') as f:
            call('qstat -x -tJ ' + jobid + ' | '
                + 'tail -n 1 ' + ' | '
                + 'awk \'{print $5}\'',
                stdout=f)

        with open(PATH.SYSTEM+'/'+'job_status', 'r') as f:
            line = f.readline()
            state = line.strip()

        return state


    ### utility function

    def _timestamp(self):
        with open(PATH.SYSTEM+'/'+'timestamps', 'a') as f:
            line = time.strftime('%H:%M:%S')+'\n'
            f.write(line)

    def save_kwargs(self, classname, method, kwargs):
        kwargspath = join(PATH.OUTPUT, 'kwargs')
        kwargsfile = join(kwargspath, classname+'_'+method+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)


