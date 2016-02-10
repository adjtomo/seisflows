
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
                + '-a intelmpi '
                + '-J %s ' % PAR.SUBTITLE
                + '-q LAURE_USERS '
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.log')
                + '-n %d ' % 16
                + '-e %s ' % (PATH.SUBMIT+'/'+'error.log')
                + '-R "span[ptile=%d' % PAR.NODESIZE + ']" '
                + '-W %d:00 ' % PAR.WALLTIME
                +  findpath('system') +'/'+ 'lsf/wrapper_bsub '
                + PATH.OUTPUT)



    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts.
        """
        self.save_objects()
        self.save_kwargs(classname, funcname, kwargs)
        jobs = self.launch(classname, funcname, hosts)
        while 1:
            # pauses python program
            time.sleep(60*PAR.SLEEPTIME)
            self.timestamp()
            isdone, jobs = self.task_status(classname, funcname, jobs)
            if isdone:
                return



    def launch(self, classname, funcname, hosts='all'):
        unix.mkdir(PATH.SYSTEM)

        # prepare bsub arguments
        if hosts == 'all':
            args = ('[%d-%d]' % (1,PAR.NSRC) + '%' + '%d" ' %(PAR.NTASK)                             
                   + '-o %s ' % (PATH.SUBMIT+'/'+'output.lsf/'+'%J_%I'))

        elif hosts == 'head':
            args = ('[%d-%d]" ' % (1,1)
                   + '-o %s ' % (PATH.SUBMIT+'/'+'output.lsf/'+'%J'))

        # submit job
        with open(PATH.SYSTEM+'/'+'job_id', 'w') as f:
            subprocess.call('bsub '
                + '-a intelmpi '
                + '-q LAURE_USERS '
                + '-n %d ' % PAR.NPROC 
                + '-R "span[ptile=%d' % PAR.NODESIZE + ']" '
                + '-W %d:00 ' % PAR.STEPTIME
                + '-J "%s' %PAR.SUBTITLE
                + args
                + findpath('system') +'/'+ 'lsf/wrapper_srun '
                + PATH.OUTPUT + ' '
                + classname + ' '
                + funcname + ' ',
                shell=1,
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
                #print msg.TaskError_SLURM % (classname, funcname, job)
                sys.exit(-1)

	# return True if all elements of the states are true (or if the state is empty). 
        isdone = all(states)

        return isdone, jobs




    def mpiargs(self):
        #return 'mpirun '
        #return ('/apps/lsf/cluster_ICEX/8.3/linux2.6-glibc2.3-x86_64/bin/mpirun.lsf '
        #        + '-genv I_MPI_EXTRA_FILESYSTEM 1 -genv I_MPI_EXTRA_FILESYSTEM_LIST lustre '
        #        + '-genv I_MPI_PIN 0 -genv I_MPI_FALLBACK 0 -_MSG_SIZE 4194304 -pam ')
        return ('/apps/lsf/cluster_ICEX/8.3/linux2.6-glibc2.3-x86_64/bin/mpirun.lsf '
		+ '-genv I_MPI_EXTRA_FILESYSTEM 1 -genv I_MPI_EXTRA_FILESYSTEM_LIST lustre '
		+ '-genv I_MPI_PIN 0 -genv I_MPI_FALLBACK 0 -genv I_MPI_RDMA_RNDV_WRITE 1 -genv I_MPI_RDMA_MAX_MSG_SIZE 4194304 -pam '
		+ ' "-n %s " ' % PAR.NPROC )


    def getstate(self, jobid):
        """ Retrives job state from LSF database
        """
        # sacct - report job (or job step) accounting info about active or completed jobs 
        with open(PATH.SYSTEM+'/'+'job_status', 'w') as f:
            # subprocess.call('sacct -n -o state -j '+jobid, shell=True, stdout=f)
            # SLURM: 
            # -n, --noheader: No heading will be added to the output
            # -o, --format: Comma separated list of fields.
            # LSF:
            # -N host_name
            #subprocess.call('bjobs -noheader -a -d ' + jobid, shell=True, stdout=f)
            subprocess.call('bjobs -a -d "' + jobid + '"', shell=True, stdout=f)
        with open(PATH.SYSTEM+'/'+'job_status', 'r') as f:
            line_buf = f.readline()
            line = f.readline()
	    state = line.split()[2].strip()
            #state = line.strip()

        return state



#    def getnode(self):
#        """ Gets number of running task
#        """
#        try:
#            return int(os.getenv('SEISFLOWS_TASK_ID'))
#        except:
#            try:
#                return int(os.getenv('SLURM_ARRAY_TASK_ID'))
#            except:
#                raise Exception("TASK_ID environment variable not defined.")

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

