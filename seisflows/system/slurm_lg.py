#
# This is Seisflows
#
# See LICENCE file
#
###############################################################################

# Import system modules
import os
import sys
import time
from os.path import abspath, basename, join
from subprocess import check_output

# Import utilitaries
import math

# Local imports
from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.tools import call, findpath, saveobj, timestamp
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class slurm_lg(custom_import('system', 'base')):
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

        # time allocated for workflow in minutes
        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        # time allocated for each individual task in minutes
        if 'TASKTIME' not in PAR:
            setattr(PAR, 'TASKTIME', 15.)

        # number of tasks
        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        # number of cores per task
        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        # limit on number of concurrent tasks
        if 'NTASKMAX' not in PAR:
            setattr(PAR, 'NTASKMAX', 100)

        # number of cores per node
        if 'NODESIZE' not in PAR:
            raise ParameterError(PAR, 'NODESIZE')

        # optional additional SLURM arguments
        if 'SLURMARGS' not in PAR:
            setattr(PAR, 'SLURMARGS', '')

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
        unix.mkdir(PATH.WORKDIR+'/'+'output.slurm')

        workflow.checkpoint()

        # prepare sbatch arguments
        call('sbatch '
             + '%s ' % PAR.SLURMARGS
             + '--job-name=%s ' % PAR.TITLE
             + '--output %s ' % (PATH.WORKDIR+'/'+'output.log')
             + '--ntasks-per-node=%d ' % PAR.NODESIZE
             + '--nodes=%d ' % 1
             + '--time=%d ' % PAR.WALLTIME
             + findpath('seisflows.system') + '/' + 'wrappers/submit '
             + PATH.OUTPUT)

    def run(self, classname, method, *args, **kwargs):
        """ Runs task multiple times in embarrassingly parallel fasion

          Executes classname.method(*args, **kwargs) NTASK times, each time on
          NPROC cpu cores
        """
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        # submit job array
        stdout = check_output(
                   'sbatch %s ' % PAR.SLURMARGS
                   + '--job-name=%s ' % PAR.TITLE
                   + '--nodes=%d ' % math.ceil(PAR.NPROC/float(PAR.NODESIZE))
                   + '--ntasks-per-node=%d ' % PAR.NODESIZE
                   + '--ntasks=%d ' % PAR.NPROC
                   + '--time=%d ' % PAR.TASKTIME
                   + '--array=%d-%d ' % (0, (PAR.NTASK-1) % PAR.NTASKMAX)
                   + '--output %s ' % (PATH.WORKDIR + '/' + 'output.slurm/' +
                                       '%A_%a')
                   + '%s ' % (findpath('seisflows.system') + '/' +
                              'wrappers/run')
                   + '%s ' % PATH.OUTPUT
                   + '%s ' % classname
                   + '%s ' % method
                   + '%s ' % PAR.ENVIRONS,
                   shell=True)

        # keep track of job ids
        jobs = self.job_id_list(stdout, PAR.NTASK)

        # check job array completion status
        while True:
            # wait a few seconds between queries
            time.sleep(5)

            isdone, jobs = self.job_array_status(classname, method, jobs)
            if isdone:
                return

    def run_single(self, classname, method, *args, **kwargs):
        """ Runs task a single time

          Executes classname.method(*args, **kwargs) a single time on NPROC
          cpu cores
        """
        self.checkpoint(PATH.OUTPUT, classname, method, args, kwargs)

        # submit job
        stdout = check_output(
                   'sbatch %s ' % PAR.SLURMARGS
                   + '--job-name=%s ' % PAR.TITLE
                   + '--nodes=%d ' % math.ceil(PAR.NPROC/float(PAR.NODESIZE))
                   + '--ntasks-per-node=%d ' % PAR.NODESIZE
                   + '--ntasks=%d ' % PAR.NPROC
                   + '--time=%d ' % PAR.TASKTIME
                   + '--array=%d-%d ' % (0, 0)
                   + '--output %s ' % (PATH.WORKDIR + '/' + 'output.slurm/' +
                                       '%A_%a')
                   + '%s ' % (findpath('seisflows.system') + '/' +
                              'wrappers/run')
                   + '%s ' % PATH.OUTPUT
                   + '%s ' % classname
                   + '%s ' % method
                   + '%s ' % PAR.ENVIRONS
                   + '%s ' % 'SEISFLOWS_TASKID=0',
                   shell=True)

        # keep track of job ids
        jobs = self.job_id_list(stdout, 1)

        # check job completion status
        while True:
            # wait a few seconds between queries
            time.sleep(5)

            isdone, jobs = self.job_array_status(classname, method, jobs)
            if isdone:
                return

    def mpiexec(self):
        """ Specifies MPI executable used to invoke solver
        """
        return 'srun '

    def taskid(self):
        """ Provides a unique identifier for each running task
        """
        try:
            return int(os.getenv('SEISFLOWS_TASKID'))
        except:
            return int(os.getenv('SLURM_ARRAY_TASK_ID'))

    # Job array methods

    def job_array_status(self, classname, method, jobs):
        """ Determines completion status of job or job array
        """
        states = []
        for job in jobs:
            state = self.job_status(job)
            if state in ['TIMEOUT']:
                print msg.TimoutError % (classname, method, job, PAR.TASKTIME)
                sys.exit(-1)
            elif state in ['FAILED', 'NODE_FAIL']:
                print msg.TaskError_SLURM % (classname, method, job)
                sys.exit(-1)
            elif state in ['COMPLETED']:
                states += [1]
            else:
                states += [0]

        isdone = all(states)

        return isdone, jobs

    def job_id_list(self, stdout, ntask):
        """ Parses job id list from sbatch standard output
        """
        job_id = stdout.split()[-1].strip()
        return [job_id+'_'+str(ii) for ii in range(ntask)]

    def job_status(self, job):
        """ Queries completion status of a single job
        """
        stdout = check_output(
            'sacct -n -o jobid,state -j ' + job.split('_')[0],
            shell=True)

        state = ''
        lines = stdout.strip().split('\n')
        for line in lines:
            if line.split()[0] == job:
                state = line.split()[1]
        return state
