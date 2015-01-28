
import os
import math
import sys
import subprocess
import time
from os.path import abspath, join

from seisflows.tools import unix
from seisflows.tools.code import saveobj
from seisflows.tools.config import findpath, ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

save_objects = OBJ.save
save_parameters = PAR.save
save_paths = PATH.save


class slurm_lg_job(object):
    """ System interface class

      Provides an interface through which to submit jobs, run tasks in serial
      or parallel, and perform other system functions.

      One of several system interface classes that together provide a consistent
      interface across different computer environemnts. Each class implements a
      standard sets of methods, hiding the details associated with, for example,
      a particular filesystem or job scheduler.
    """


    def check(self):
        """ Checks parameters and paths
        """

        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', unix.basename(abspath('.')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR, 'SUBTITLE', unix.basename(abspath('..')))

        # check parameters
        if 'NTASK' not in PAR:
            raise Exception

        if 'NPROC' not in PAR:
            raise Exception

        if 'NPROC_PER_NODE' not in PAR:
            raise Exception

        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        if 'STEPTIME' not in PAR:
            setattr(PAR, 'STEPTIME', 30.)

        if 'SLEEPTIME' not in PAR:
            PAR.SLEEPTIME = 1.

        if 'RETRY' not in PAR:
            PAR.RETRY = False

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        # check paths
        if 'GLOBAL' not in PATH:
            setattr(PATH, 'GLOBAL', join(abspath('.'), 'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', '')

        if 'SUBMIT' not in PATH:
            setattr(PATH, 'SUBMIT', unix.pwd())

        if 'OUTPUT' not in PATH:
            setattr(PATH, 'OUTPUT', join(PATH.SUBMIT, 'output'))

        if 'SYSTEM' not in PATH:
            setattr(PATH, 'SYSTEM', join(PATH.GLOBAL, 'system'))


    def submit(self, workflow):
        """ Submits job
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)

        # save current state
        save_objects('SeisflowsObjects')
        save_parameters('SeisflowsParameters.json')
        save_paths('SeisflowsPaths.json')

        unix.mkdir(PATH.SUBMIT+'/'+'output.slurm')

        args = ('sbatch '
          + '--job-name=%s ' % PAR.TITLE
          + '--output %s ' % (PATH.SUBMIT+'/'+'output.log')
          + '--ntasks-per-node=%d ' % PAR.NPROC_PER_NODE
          + '--nodes=%d ' % 1
          + '--time=%d ' % PAR.WALLTIME
          + findpath('system') +'/'+ 'slurm/wrapper_sbatch '
          + PATH.OUTPUT)

        subprocess.call(args, shell=1)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts
        """
        unix.mkdir(PATH.SYSTEM)

        if PAR.VERBOSE >= 2:
            print 'running', funcname

        # save current state
        save_objects(join(PATH.OUTPUT, 'SeisflowsObjects'))

        # save keyword arguments
        kwargspath = join(PATH.OUTPUT, 'SeisflowsObjects', classname+'_kwargs')
        kwargsfile = join(kwargspath, funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)

        if hosts == 'all':
            # run on all available nodes
            jobs = self.launch(classname, funcname, hosts)
            while 1:
                time.sleep(60.*PAR.SLEEPTIME)
                self.timestamp()
                isdone, jobs = self.task_status(classname, funcname, jobs)
                if isdone:
                    return
        elif hosts == 'head':
            # run on head node
            jobs = self.launch(classname, funcname, hosts)
            while 1:
                time.sleep(60.*PAR.SLEEPTIME)
                self.timestamp()
                isdone, jobs = self.task_status(classname, funcname, jobs)
                if isdone:
                    return
        else:
            raise Exception


    def launch(self, classname, funcname, hosts='all'):
        # prepare sbatch arguments
        sbatch = 'sbatch '                                                   \
            + '--job-name=%s ' % PAR.TITLE                                   \
            + '--nodes=%d ' % math.ceil(PAR.NPROC/float(PAR.NPROC_PER_NODE)) \
            + '--ntasks-per-node=%d ' % PAR.NPROC_PER_NODE                   \
            + '--time=%d ' % PAR.STEPTIME

        args = findpath('system') +'/'+ 'slurm/wrapper_srun '        \
            + PATH.OUTPUT + ' '                                      \
            + classname + ' '                                        \
            + funcname + ' '

        # call sbatch
        if hosts == 'all':
            with open(PATH.SYSTEM+'/'+'job_id', 'w') as f:
                subprocess.call(sbatch
                  + '--array=%d-%d ' % (0,PAR.NTASK-1)
                  + '--output %s ' % (PATH.SUBMIT+'/'+'output.slurm/'+'%A_%a')
                  + args,
                  shell=1, stdout=f)
        elif hosts == 'head':
            with open(PATH.SYSTEM+'/'+'job_id', 'w') as f:
                subprocess.call(sbatch
                  + '--array=%d-%d ' % (0,0)
                  + '--output %s ' % (PATH.SUBMIT+'/'+'output.slurm/'+'%A_%a')
                  + args,
                  shell=1, stdout=f)

        # collect job info
        with open(PATH.SYSTEM+'/'+'job_id', 'r') as f:
            line = f.readline()

        job_id = line.split()[-1].strip()
        job_ids = [job_id+'_'+str(id) for id in range(self.mylen(hosts))]
        task_ids = range(self.mylen(hosts))

        return dict(zip(job_ids, task_ids))


    def task_status(self, classname, funcname, jobs):
        job_ids = jobs.keys()
        task_ids = jobs.values()

        # query job database
        with open(PATH.SYSTEM+'/'+'job_status','w') as f:
            subprocess.call('sacct -n -o jobid,state', shell=1, stdout=f)

        with open(PATH.SYSTEM+'/'+'job_status', 'r') as f:
            lines = f.readlines()

        isdone = 1
        for line in lines:
            job_id, status = line.split()
            if job_id not in job_ids:
                continue
            # check job status
            task_id = jobs[job_id]
            if status in ['PENDING', 'RUNNING']:
                isdone = 0
            elif status in ['FAILED']:
                print 'JOB FAILED:', job_id
                if PAR.RETRY:
                    jobs.__delitem__(job_id)
                    job = self.launch(classname, funcname, hosts=task_id)
                    job_id = job.keys().pop()
                    jobs.__setitem__(job_id, task_id)
                    isdone = 0
                    print 'RELAUNCHING...'
                    print ''
                else:
                    raise Exception
        return isdone, jobs


    def mpiargs(self):
        return 'srun '

    def getnode(self):
        """ Gets number of running task
        """
        try:
            return int(os.getenv('SLURM_ARRAY_TASK_ID'))
        except:
            raise Exception
    def mylen(self, hosts):
        if hosts=='all':
            return PAR.NTASK
        else:
            return 1

    def timestamp(self):
        with open(PATH.SYSTEM+'/'+'timestamps', 'a') as f:
            line = time.strftime('%H:%M:%S')+'\n'
            f.write(line)

