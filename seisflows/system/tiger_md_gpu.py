
import os
import sys
from os.path import abspath, join

from seisflows.tools import unix
from seisflows.tools.code import call, findpath, saveobj
from seisflows.tools.config import ParameterError, custom_import, \
    SeisflowsObjects, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class tiger_md_gpu(custom_import('system', 'tiger_md')):
    """ An interface through which to submit workflows, run tasks in serial or 
      parallel, and perform other system functions.

      By hiding environment details behind a python interface layer, these 
      classes provide a consistent command set across different computing
      environments.

      For important additional information, please see 
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration
    """

    def check(self):
        """ Checks parameters and paths
        """
        raise NotImplementedError('Provided by Etienne Bachmann. Not recently testested and not likely to work out of the box.')

        # why does Etienne have it this way?
        if 'NGPU' not in PAR:
            setattr(PAR, 'NGPU', 4)

        super(tiger_md_gpu, self).check()


    def submit(self, workflow):
        """ Submits workflow
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)

        self.checkpoint()

        if not exists(PATH.SUBMIT + '/' + 'scratch'):
            unix.ln(PATH.SCRATCH, PATH.SUBMIT + '/' + 'scratch')

        call('sbatch '
                + '--job-name=%s ' % PAR.SUBTITLE
                + '--output=%s ' % (PATH.SUBMIT +'/'+ 'output.log')
                + '--nodes 1 '
                + '--ntasks=% ' % PAR.NGPU
                + '--ntasks-per-socket=%d ' % PAR.NGPU
                + '--gres=gpu:%d ' % PAR.NGPU
                + '--time=%d ' % PAR.WALLTIME
                + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts
        """
        self.checkpoint()
        self.save_kwargs(classname, funcname, kwargs)


        if hosts == 'all':
            call('srun '
                    + '--wait=0 '
                    + join(findpath('seisflows.system'), 'wrappers/run ')
                    + PATH.OUTPUT + ' '
                    + classname + ' '
                    + funcname)

        elif hosts == 'head':
            # run on head node
            call('srun '
                    + '--wait=0 '
                    + join(findpath('seisflows.system'), 'wrappers/run_head ')
                    + PATH.OUTPUT + ' '
                    + classname + ' '
                    + funcname)

    def getnode(self):
        """ Gets number of running task
        """
        gid = os.getenv('SLURM_GTIDS').split(',')
        lid = int(os.getenv('SLURM_LOCALID'))
        return int(gid[lid])


    def mpiexec(self):
        return 'mpirun -np %d '%PAR.NPROC


    def save_kwargs(self, classname, funcname, kwargs):
        kwargspath = join(PATH.OUTPUT, 'SeisflowsObjects', classname+'_kwargs')
        kwargsfile = join(kwargspath, funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)

