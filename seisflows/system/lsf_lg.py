
import os
import math
import sys
import subprocess
import time
from os.path import abspath, join

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.code import saveobj
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, findpath, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class lsf_lg(loadclass('system', 'base')):
    """ An interface through which to submit workflows, run tasks in serial or 
      parallel, and perform other system functions.

      By hiding environment details behind a python interface layer, these 
      classes provide a consistent command set across different computing
      environments.

      For more informations, see 
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-interfaces
    """

    def check(self):
        """ Checks parameters and paths
        """

        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', unix.basename(abspath('..')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR, 'SUBTITLE', unix.basename(abspath('.')))

        # check parameters
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

        if 'LSF_ARGS' not in PAR:
            setattr(PAR, 'LSF_ARGS', '')

        # check paths
        if 'GLOBAL' not in PATH:
            setattr(PATH, 'GLOBAL', join(abspath('.'), 'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SUBMIT' not in PATH:
            setattr(PATH, 'SUBMIT', unix.pwd())

        if 'OUTPUT' not in PATH:
            setattr(PATH, 'OUTPUT', join(PATH.SUBMIT, 'output'))

        if 'SYSTEM' not in PATH:
            setattr(PATH, 'SYSTEM', join(PATH.GLOBAL, 'system'))


    def submit(self, workflow):
        """ Submits workflow
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)
        unix.mkdir(PATH.SUBMIT+'/'+'output.lsf')

        self.save_objects()
        self.save_parameters()
        self.save_paths()

        # prepare bsub arguments
        unix.run('bsub '
                + PAR.LSF_ARGS + ' '
                + '-J %s ' % PAR.SUBTITLE
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.log')
                + '-n %d ' % PAR.NODESIZE
                + '-e %s ' % (PATH.SUBMIT+'/'+'error.log')
                + '-R "span[ptile=%d' % PAR.NODESIZE + ']" '
                + '-W %d:00 ' % PAR.WALLTIME
                +  findpath('system') +'/'+ 'wrappers/submit '
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
            subprocess.call('bsub '
                + PAR.LSF_ARGS + ' '
                + '-n %d ' % PAR.NPROC 
                + '-R "span[ptile=%d]" ' % PAR.NODESIZE
                + '-W %d:00 ' % PAR.STEPTIME
                + '-J "%s' %PAR.SUBTITLE
                + self.launch_args(hosts)
                + findpath('system') +'/'+ 'wrapper/run '
                + PATH.OUTPUT + ' '
                + classname + ' '
                + funcname + ' ',
                shell=True,
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
            args += '-o %s ' % (PATH.SUBMIT+'/'+'output.lsf/'+'%J_%I'))

        elif hosts == 'head':
            args = ''
            args += '[%d-%d]' % (1, 1)
            args += '-o %s ' % (PATH.SUBMIT+'/'+'output.lsf/'+'%J'))

        return args



    def mpiargs(self):
        return 'mpirun '


    def getstate(self, jobid):
        """ Retrives job state from LSF database
        """
        with open(PATH.SYSTEM+'/'+'job_status', 'w') as f:
            subprocess.call('bjobs -a -d "' + jobid + '"', shell=True, stdout=f)
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

