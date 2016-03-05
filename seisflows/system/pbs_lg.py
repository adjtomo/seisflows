
import os
import math
import sys
import subprocess
import time
from os.path import abspath, join

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.code import findpath, saveobj
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class pbs_lg(loadclass('system', 'base')):
    """ An interface through which to submit workflows, run tasks in serial or
      parallel, and perform other system functions.

      By hiding environment details behind a python interface layer, these
      classes provide a consistent command set across different computing
      environments.

      Intermediate files are written to a global scratch path PATH.SCRATCH,
      which must be accessible to all compute nodes.

      Optionally, users can provide a local scratch path PATH.LOCAL if each
      compute node has its own local filesystem.

      For more informations, see
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-interfaces
    """

    def check(self):
        """ Checks parameters and paths
        """

        # check parameters
        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', unix.basename(abspath('.')))

        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        if 'STEPTIME' not in PAR:
            setattr(PAR, 'STEPTIME', 30.)

        if 'SLEEPTIME' not in PAR:
            setattr(PAR, 'SLEEPTIME', 1.)

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        if 'NODESIZE' not in PAR:
            raise ParameterError(PAR, 'NODESIZE')

        if 'PBS_ARGS' not in PAR:
            setattr(PAR, 'PBS_ARGS', '')

        # check paths
        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH', join(abspath('.'), 'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SUBMIT' not in PATH:
            setattr(PATH, 'SUBMIT', unix.pwd())

        if 'OUTPUT' not in PATH:
            setattr(PATH, 'OUTPUT', join(PATH.SUBMIT, 'output'))

        if 'SYSTEM' not in PATH:
            setattr(PATH, 'SYSTEM', join(PATH.SCRATCH, 'system'))


    def submit(self, workflow):
        """ Submits workflow
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)
        unix.mkdir(PATH.SUBMIT+'/'+'output.pbs')

        self.checkpoint()

        hours = PAR.WALLTIME/60
        minutes = PAR.WALLTIME%60
        walltime = 'walltime=%02d:%02d:00 ' % (hours, minutes)

        ncpus = PAR.NODESIZE
        mpiprocs = PAR.NODESIZE

        # prepare qsub arguments
        unix.run( 'qsub '
                + PAR.PBS_ARGS + ' '
                + '-l select=1:ncpus=%d:mpiprocs=%d ' % (ncpus, mpiprocs)
                + '-l %s ' % walltime
                + '-N %s ' % PAR.TITLE
                + '-j %s '%'oe'
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.log')
                + '-V '
                + ' -- ' + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """ Runs tasks in serial or parallel on specified hosts.
        """
        self.checkpoint()

        self.save_kwargs(classname, funcname, kwargs)
        jobs = self._launch(classname, funcname, hosts)
        while True:
            time.sleep(60.*PAR.SLEEPTIME)
            self._timestamp()
            isdone, jobs = self._status(classname, funcname, jobs)
            if isdone:
                return


    def mpiargs(self):
        return 'mpirun '


    def getnode(self):
        """ Gets number of running task
        """
        try:
            return os.getenv('PBS_ARRAY_INDEX')
        except:
            raise Exception("PBS_ARRAY_INDEX environment variable not defined.")


    ### private methods

    def _launch(self, classname, funcname, hosts='all'):
        unix.mkdir(PATH.SYSTEM)

        nodes = math.ceil(PAR.NTASK/float(PAR.NODESIZE))
        ncpus = PAR.NPROC
        mpiprocs = PAR.NPROC

        hours = PAR.STEPTIME/60
        minutes = PAR.STEPTIME%60
        walltime = 'walltime=%02d:%02d:00 '%(hours, minutes)

        # submit job
        with open(PATH.SYSTEM+'/'+'job_id', 'w') as f:
            subprocess.call('qsub '
                + PAR.PBS_ARGS + ' '
                + '-l select=%d:ncpus=%d:mpiprocs=%d ' (nodes, ncpus, mpiprocs)
                + '-l %s ' % walltime
                + '-J 0-%s ' % (PAR.NTASK-1)
                + '-N %s ' % PAR.TITLE
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.pbs/' + '$PBS_ARRAYID')
                + '-r y '
                + '-j oe '
                + '-V '
                + self.launch_args(hosts)
                + PATH.OUTPUT + ' '
                + classname + ' '
                + funcname + ' '
                + findpath('seisflows'),
                stdout=f,
                shell=True)
           
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


    def launch_args(self, hosts):
        if hosts == 'all':
          arg = ('-J 0-%s ' % (PAR.NTASK-1)
                +'-o %s ' % (PATH.SUBMIT+'/'+'output.pbs/' + '$PBS_ARRAYID')
                + ' -- ' + findpath('seisflows.system') +'/'+ 'wrappers/run_pbsdsh ')

        elif hosts == 'head':
          arg = ('-J 0-0 '
                 +'-o %s ' % (PATH.SUBMIT+'/'+'output.pbs/' + '$PBS_JOBID')
                 + ' -- ' + findpath('seisflows.system') +'/'+ 'wrappers/run_pbsdsh_head ')

        return arg


    def _status(self, classname, funcname, jobs):
        """ Determines completion status of one or more jobs
        """
        for job in jobs:
            state = self._query(job)
            states = []
            if state in ['C']:
                states += [1]
            else:
                states += [0]
            if state in ['F']:
                print msg.TaskError_PBS % (classname, funcname, job)
                sys.exit(-1)
        isdone = all(states)

        return isdone, jobs


    def _query(self, jobid):
        """ Queries job state from PBS database
        """
        # TODO: replace shell utilities with native Python
        with open(PATH.SYSTEM+'/'+'job_status', 'w') as f:
            subprocess.call('qstat -x -tJ ' + jobid + ' | '
                + 'tail -n 1 ' + ' | '
                + 'awk \'{print $5}\'',
                shell=True,
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

    def save_kwargs(self, classname, funcname, kwargs):
        kwargspath = join(PATH.OUTPUT, 'SeisflowsObjects', classname+'_kwargs')
        kwargsfile = join(kwargspath, funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)


