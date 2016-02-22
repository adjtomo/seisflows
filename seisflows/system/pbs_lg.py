
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


class pbs_lg(loadclass('system', 'base')):
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
        unix.mkdir(PATH.SUBMIT+'/'+'output.pbs')

        self.checkpoint()

        # construct resource list
        nodes = PAR.NTASK/PAR.NODESIZE
        cores = PAR.NTASK%PAR.NODESIZE
        hours = PAR.WALLTIME/60
        minutes = PAR.WALLTIME%60
        resources = 'walltime=%02d:%02d:00 '%(hours, 10)
        #if nodes == 0:
        #    resources += ',nodes=1:ppn=%d'%cores
        #lif cores == 0:
         #   resources += ',nodes=%d:ppn=16'%nodes
        #else:
        #    resources += ',nodes=%d:ppn=16+1:ppn=%d'%(nodes, cores)

        # prepare qsub arguments
        print('qsub -l select=1:ncpus=32:mpiprocs=2 -A ERDCH38424KSC -q debug -N %s -o %s -l %s -j %s')%(PAR.SUBTITLE, 'output.log', resources, 'oe'+findpath('system')+'/'+'wrapper/submit ' + PATH.OUTPUT)

        args = ('qsub '
                + '-l select=1:ncpus=32:mpiprocs=4 '
                + '-l %s '%resources
                + '-q standard '
                + '-A ERDCH38424KSC '
                + '-N %s ' % PAR.SUBTITLE
                + '-j %s '%'oe'
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.log')
                + '-V '
                + ' -- ' + findpath('system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)

        print 'args:', args

        subprocess.call(args, shell=True)

    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts.
        """
        print("made it to beginning of run")
        self.checkpoint()

        self.save_kwargs(classname, funcname, kwargs)
        jobs = self._launch(classname, funcname, hosts)
        while 1:
            time.sleep(10.*PAR.SLEEPTIME)
            self._timestamp()
            isdone, jobs = self._status(classname, funcname, jobs)
            if isdone:
                return


    def mpiargs(self):
        return 'aprun -n %d' % 1


    def getnode(self):
        """ Gets number of running task
        """
        try:
            return os.getenv('PBS_ARRAY_INDEX')
        except:
            raise Exception("TASK_ID environment variable not defined.")


    ### private methods

    def _launch(self, classname, funcname, hosts='all'):
        unix.mkdir(PATH.SYSTEM)

        # submit job
        with open(PATH.SYSTEM+'/'+'job_id', 'w') as f:
            hours = PAR.WALLTIME/60
            minutes = PAR.WALLTIME%60
            resources = 'walltime=%02d:%02d:00 '%(hours, minutes)
            args = ('/opt/pbs/12.1.1.131502/bin/qsub '
                + '-l select=1:ncpus=32:mpiprocs=4 '
                + '-l %s '%resources
                + '-q standard '
                + '-A ERDCH38424KSC '
                + '-J 0-%s ' % (PAR.NTASK-1)
                + '-N %s ' % PAR.TITLE
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.pbs/' + '$PBS_ARRAYID')
                + '-r y '
                + '-j oe '
                + '-V '
                + ' -- ' + findpath('system') +'/'+ 'wrappers/run_pbsdsh '
                + PATH.OUTPUT + ' '
                + classname + ' '
                + funcname + ' '
                + '/work/jas11/SEISFLOWS')
           
            print(args)
            subprocess.call(args, shell=1, stdout=f)

        print("made it here to job id")
        # retrieve job ids
        with open(PATH.SYSTEM+'/'+'job_id', 'r') as f:
            line = f.readline()
            job = line.split()[-1].strip()
        if hosts == 'all' and PAR.NTASK > 1:
            nn = range(PAR.NTASK)
            # take number[].sdb and replace with number[str(ii)]].sdb
            jobMain = job.split('[',1)[0]
            print(jobMain)
            return [jobMain+'['+str(ii)+'].sdb' for ii in nn]
        else:
            return [job]

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
                #print msg.TaskError_SLURM % (classname, funcname, job)
                sys.exit(-1)
        isdone = all(states)

        return isdone, jobs


    def _query(self, jobid):
        """ Queries job state from PBS database
        """
        with open(PATH.SYSTEM+'/'+'job_status', 'w') as f:
            subprocess.call('/opt/pbs/12.1.1.131502/bin/qstat -x -tJ ' +jobid+ ' | tail -n 1 | awk \'{print $5}\'', shell=True, stdout=f)

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


