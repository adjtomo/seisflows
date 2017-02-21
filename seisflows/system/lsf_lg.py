
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
        # name of job
        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', basename(abspath('.')))

        # time allocated for entire workflow
        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        # time allocated for each individual task
        if 'STEPTIME' not in PAR:
            setattr(PAR, 'STEPTIME', 15.)

        # number of tasks
        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        # number of cores per task
        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

         # number of cores per node
        if 'NODESIZE' not in PAR:
            raise ParameterError(PAR, 'NODESIZE')

        # optional additional LSF arguments
        if 'LSFARGS' not in PAR:
            setattr(PAR, 'LSFARGS', '')

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
        unix.mkdir(PATH.WORKDIR+'/'+'output.lsf')

        self.checkpoint()

        # prepare bsub arguments
        call('bsub '
                + '%s ' % PAR.LSFARGS
                + '-J %s ' % PAR.TITLE
                + '-o %s ' % (PATH.WORKDIR+'/'+'output.log')
                + '-n %d ' % PAR.NODESIZE
                + '-e %s ' % (PATH.WORKDIR+'/'+'error.log')
                + '-R "span[ptile=%d]" ' % PAR.NODESIZE
                + '-W %d:00 ' % PAR.WALLTIME
                +  findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts.
        """
        self.save_objects()
        self.save_kwargs(classname, funcname, kwargs)
        jobs = self.submit_job_array(classname, funcname, hosts)
        while True:
            # wait 30 seconds before checking status again
            time.sleep(30)

            self.timestamp()
            isdone, jobs = self.job_status(classname, funcname, jobs)
            if isdone:
                return


    def submit_job_array(self, classname, funcname, hosts='all'):
        # submit job
        with open(PATH.SYSTEM+'/'+'job_id', 'w') as f:
            call(self.job_array_cmd(classname, funcname, hosts), stdout=f)

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


    def job_array_cmd(self, classname, funcname, hosts):
        return ('bsub '
            + '%s ' % PAR.LSFARGS
            + '-n %d ' % PAR.NPROC
            + '-R "span[ptile=%d]" ' % PAR.NODESIZE
            + '-W %d:00 ' % PAR.STEPTIME
            + '-J "%s' %PAR.TITLE
            + self.launch_args(hosts)
            + findpath('seisflows.system') +'/'+ 'wrapper/run '
            + PATH.OUTPUT + ' '
            + classname + ' '
            + funcname + ' '
            + PAR.ENVIRONS)


    def job_array_args(self, hosts):
        if hosts == 'all':
            args = ''
            args += '[%d-%d] %% %d' % (1, PAR.NSRC, PAR.NTASK)
            args += '-o %s ' % (PATH.WORKDIR+'/'+'output.lsf/'+'%J_%I')

        elif hosts == 'head':
            args = ''
            args += '[%d-%d]' % (1, 1)
            args += '-o %s ' % (PATH.WORKDIR+'/'+'output.lsf/'+'%J')

        return args


    def job_status(self, classname, funcname, jobs):
        # query lsf database
        for job in jobs:
            state = self._query(job)
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


    def mpiexec(self):
        """ Specifies MPI exectuable; used to invoke solver
        """
        return 'mpiexec '


    def _query(self, jobid):
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
        kwargspath = join(PATH.OUTPUT, 'kwargs')
        kwargsfile = join(kwargspath, classname+'_'+funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)

