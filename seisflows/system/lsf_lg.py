
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


class lsf_lg(custom_import('system', 'base')):
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

        # check parameters
        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', basename(abspath('.')))

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

        if 'LSFARGS' not in PAR:
            setattr(PAR, 'LSFARGS', '')

        # check paths
        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH', join(abspath('.'), 'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SUBMIT' not in PATH:
            setattr(PATH, 'SUBMIT', abspath('.'))

        if 'OUTPUT' not in PATH:
            setattr(PATH, 'OUTPUT', join(PATH.SUBMIT, 'output'))

        if 'SYSTEM' not in PATH:
            setattr(PATH, 'SYSTEM', join(PATH.SCRATCH, 'system'))


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
                + '-n %d ' % PAR.NODESIZE
                + '-e %s ' % (PATH.SUBMIT+'/'+'error.log')
                + '-R "span[ptile=%d]" ' % PAR.NODESIZE
                + '-W %d:00 ' % PAR.WALLTIME
                +  findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts.
        """
        self.save_objects()
        self.save_kwargs(classname, funcname, kwargs)
        jobs = self.launch(classname, funcname, hosts)
        while True:
            # wait a few seconds before checking status
            time.sleep(60*PAR.SLEEPTIME)

            self.timestamp()
            isdone, jobs = self.task_status(classname, funcname, jobs)
            if isdone:
                return


    def launch(self, classname, funcname, hosts='all'):
        unix.mkdir(PATH.SYSTEM)

        # submit job
        with open(PATH.SYSTEM+'/'+'job_id', 'w') as f:
            call('bsub '
                + '%s ' % PAR.LSFARGS
                + '-n %d ' % PAR.NPROC 
                + '-R "span[ptile=%d]" ' % PAR.NODESIZE
                + '-W %d:00 ' % PAR.STEPTIME
                + '-J "%s' %PAR.TITLE
                + self.launch_args(hosts)
                + findpath('seisflows.system') +'/'+ 'wrapper/run '
                + PATH.OUTPUT + ' '
                + classname + ' '
                + funcname + ' ',
                stdout=f)

        # retrieve job ids
        with open(PATH.SYSTEM+'/'+'job_id', 'r') as f:
            # reads one entire line from the file
            line = f.readline()
            job_buf = line.split()[1].strip()
            job = job_buf[1:-1]
        if hosts == 'all' and PAR.NSRC > 1:
            nn = range(1,PAR.NSRC+1)
            #return [job+'_'+str(ii) for ii in nn]
            return [job+'['+str(ii)+']' for ii in nn]
        else:
            return [job]


    def task_status(self, classname, funcname, jobs):
        # query lsf database
        for job in jobs:
            state = self.getstate(job)
            states = []
            if state in ['DONE']:
                states += [1]
            else:
                states += [0]
            if state in ['EXIT']:
                print 'LSF job failed: %s ' %job
                print msg.TaskError_LSF % (classname, funcname, job)
                sys.exit(-1)
        isdone = all(states)

        return isdone, jobs


    def launch_args(self, hosts):
        if hosts == 'all':
            args = ''
            args += '[%d-%d] %% %d' % (1, PAR.NSRC, PAR.NTASK)
            args += '-o %s ' % (PATH.SUBMIT+'/'+'output.lsf/'+'%J_%I')

        elif hosts == 'head':
            args = ''
            args += '[%d-%d]' % (1, 1)
            args += '-o %s ' % (PATH.SUBMIT+'/'+'output.lsf/'+'%J')

        return args



    def mpiexec(self):
        """ Specifies MPI exectuable; used to invoke solver
        """
        return 'mpiexec '


    def getstate(self, jobid):
        """ Retrives job state from LSF database
        """
        with open(PATH.SYSTEM+'/'+'job_status', 'w') as f:
            call('bjobs -a -d "' + jobid + '"', stdout=f)
        with open(PATH.SYSTEM+'/'+'job_status', 'r') as f:
            lines = f.readlines()
            state = lines[1].split()[2].strip()
        return state


    def getnode(self):
        """ Gets number of running task
        """
        return int(os.getenv('LSB_JOBINDEX'))-1


    def timestamp(self):
        with open(PATH.SYSTEM+'/'+'timestamps', 'a') as f:
            line = time.strftime('%H:%M:%S')+'\n'
            f.write(line)


    def save_kwargs(self, classname, funcname, kwargs):
        kwargspath = join(PATH.OUTPUT, 'SeisflowsObjects', classname+'_kwargs')
        kwargsfile = join(kwargspath, funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)

