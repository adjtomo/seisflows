
import os
import math
import pickle
import socket
import sys
import subprocess
import time

from seisflows.tools import unix
from seisflows.tools.codetools import abspath, join, savejson, saveobj
from seisflows.tools.configtools import getmodule, getpath
from seisflows.tools.configtools import ParameterObject


PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')


class slurm_big_job(object):
    """ System interface class

      Provides an interface through which to submit jobs, run tasks in serial
      or parallel, and perform other system functions.

      One of several system interface classes included in SEISFLOWS to provide
      a consistent interface across different computer environemnts. Each class
      implements a standard sets of methods, hiding the details of submitting and
      running and jobs in a particular environment.
    """


    def __init__(self):
        """ Class constructor
        """

        # check parameters
        if 'NTASK' not in PAR:
            raise Exception

        if 'NPROC' not in PAR:
            raise Exception

        if 'CPUS_PER_NODE' not in PAR:
            raise Exception

        if 'WALLTIME' not in PAR:
            setattr(PAR,'WALLTIME',30.)

        if 'STEPTIME' not in PAR:
            setattr(PAR,'STEPTIME',30.)

        if 'VERBOSE' not in PAR:
            setattr(PAR,'VERBOSE',1)

        if 'TITLE' not in PAR:
            setattr(PAR,'TITLE',unix.basename(abspath('.')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR,'SUBTITLE',unix.basename(abspath('..')))

        if 'RETRY' not in PAR:
            PAR.RETRY = False

        if 'SLEEP' not in PAR:
            PAR.SLEEP = 30.

        # check paths
        if 'GLOBAL' not in PATH:
            setattr(PATH,'GLOBAL',join(abspath('.'),'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH,'LOCAL','')

        if 'SYSTEM' not in PATH:
            setattr(PATH,'SYSTEM',join(PATH.GLOBAL,'system'))

        if 'SUBMIT' not in PATH:
            setattr(PATH,'SUBMIT',unix.pwd())

        if 'SUBMIT_HOST' not in PATH:
            setattr(PATH,'SUBMIT_HOST',unix.hostname())


    def submit(self,workflow):
        """ Submits job
        """
        unix.cd(PATH.SUBMIT)

        # store parameters
        savejson(join(PATH.SUBMIT,'parameters.p'),PAR.vars)
        savejson(join(PATH.SUBMIT,'paths.p'),PATH.vars)

        subprocess.call('sbatch '
          + '--job-name=%s ' % PAR.TITLE
          + '--output %s ' % (PATH.SUBMIT+'/'+'output.log')
          + '--ntasks-per-node=%d ' % PAR.CPUS_PER_NODE
          + '--nodes=%d ' % 1
          + '--time=%d ' % PAR.WALLTIME
          + getpath('system') +'/'+ 'slurm/wrapper_sbatch '
          + PATH.SUBMIT + ' '
          + getmodule(workflow),
          shell=1)


    def run(self,task,hosts='all',**kwargs):
        """  Runs tasks in serial or parallel on specified hosts
        """
        name = task.__name__

        if PAR.VERBOSE >= 2:
            print 'running',name

        # prepare function arguments
        unix.mkdir(PATH.SYSTEM)
        file = PATH.SYSTEM+'/'+name+'.p'
        saveobj(file,kwargs)

        if hosts == 'all':
            # run on all available nodes
            jobs = self.launch(task,hosts)
            while 1:
                time.sleep(PAR.SLEEP)
                self.timestamp()
                isdone,jobs = self.task_status(task,jobs)
                if isdone:
                    return
        elif hosts == 'head':
            # run on head node
            jobs = self.launch(task,hosts)
            while 1:
                time.sleep(PAR.SLEEP)
                self.timestamp()
                isdone,jobs = self.task_status(task,jobs)
                if isdone:
                    return
        else:
            raise Exception


    def launch(self,task,hosts='all'):
        # prepare sbatch arguments
        sbatch = 'sbatch '                                           \
            + '--job-name=%s ' % PAR.TITLE                           \
            + '--output %s ' % (PATH.SYSTEM+'/'+'output.log'+'%a')   \
            + '--nodes=%d ' % math.ceil(PAR.NPROC/float(PAR.CPUS_PER_NODE)) \
            + '--ntasks-per-node=%d ' % PAR.CPUS_PER_NODE            \
            + '--time=%d ' % PAR.STEPTIME

        args = getpath('system') +'/'+ 'slurm/wrapper_srun '         \
            + PATH.SUBMIT + ' '                                  \
            + getmodule(task) + ' '                                  \
            + task.__name__ + ' '

        # call sbatch
        if hosts == 'all':
            with open(PATH.SYSTEM+'/'+'job_id','w') as f:
                subprocess.call(sbatch
                  + '--array=%d-%d ' % (0,PAR.NTASK-1)
                  + '--output %s ' % (PATH.SYSTEM+'/'+'output.log'+'%a')
                  + args,
                  shell=1,stdout=f)
        elif hosts == 'head':
            with open(PATH.SYSTEM+'/'+'job_id','w') as f:
                subprocess.call(sbatch
                  + '--export='+'MY_JOB_ID=0 '
                  + '--output %s ' % (PATH.SYSTEM+'/'+'output.log'+'0')
                  + args,
                  shell=1,stdout=f)
        else:
            itask = str(hosts)
            with open(PATH.SYSTEM+'/'+'job_id','w') as f:
                subprocess.call(sbatch
                  + '--export='+'MY_JOB_ID'+'='+itask + ' '
                  + '--output %s ' % (PATH.SYSTEM+'/'+'output.log'+itask)
                  + args,
                  shell=1,stdout=f)

        # collect job info
        with open(PATH.SYSTEM+'/'+'job_id','r') as f:
            line = f.readline()
        job_id = int(line.split()[-1].strip())

        job_ids = [str(id) for id in range(job_id,job_id+self.mylen(hosts))]
        task_ids = range(self.mylen(hosts))

        return dict(zip(job_ids,task_ids))


    def task_status(self,task,jobs):
        job_ids = jobs.keys()
        task_ids = jobs.values()

        # query job database
        with open(PATH.SYSTEM+'/'+'job_status','w') as f:
            subprocess.call('sacct -n -o jobid,state',
              shell=1,stdout=f)
        with open(PATH.SYSTEM+'/'+'job_status','r') as f:
            lines = f.readlines()

        isdone = 1
        for line in lines:
            job_id,status = line.split()
            if job_id not in job_ids:
                continue
            # check job status
            task_id = jobs[job_id]
            if status in ['PENDING','RUNNING']:
                isdone = 0
            elif status in ['FAILED']:
                print 'JOB FAILED:', job_id
                if PAR.RETRY:
                    jobs.__delitem__(job_id)
                    job = self.launch(task,hosts=task_id)
                    job_id = job.keys().pop()
                    jobs.__setitem__(job_id,task_id)
                    isdone = 0
                    print 'RELAUNCHING...'
                    print ''
                else:
                    raise Exception
        return isdone,jobs


    def mpiargs(self):
        return 'srun '


    def getnode(self):
        """ Gets number of running task
        """
        try:
            return int(os.getenv('SLURM_ARRAY_TASK_ID'))
        except:
            return int(os.getenv('MY_JOB_ID'))

    def mylen(self,hosts):
        if hosts=='all':
            return PAR.NTASK
        else:
            return 1

    def timestamp(self):
        with open(PATH.SYSTEM+'/'+'timestamps','a') as f:
            line = time.strftime('%H:%M:%S')+'\n'
            f.write(line)
