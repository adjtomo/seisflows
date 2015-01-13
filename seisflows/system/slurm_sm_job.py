
import os
import pickle
import sys
import subprocess
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


class slurm_sm_job(object):
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

        # check parameters
        if 'NTASK' not in PAR:
            raise Exception

        if 'NPROC' not in PAR:
            raise Exception

        if 'WALLTIME' not in PAR:
            setattr(PAR,'WALLTIME',30.)

        if 'VERBOSE' not in PAR:
            setattr(PAR,'VERBOSE',1)

        if 'TITLE' not in PAR:
            setattr(PAR,'TITLE',unix.basename(abspath('.')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR,'SUBTITLE',unix.basename(abspath('..')))

        # check paths
        if 'GLOBAL' not in PATH:
            setattr(PATH,'GLOBAL',join(abspath('.'),'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH,'LOCAL','')

        if 'SUBMIT' not in PATH:
            setattr(PATH,'SUBMIT',unix.pwd())

        if 'OUTPUT' not in PATH:
            setattr(PATH,'OUTPUT',join(PATH.SUBMIT,'output'))


    def submit(self, workflow):
        """ Submits job
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)

        # save current state
        save_objects('SeisflowsObjects')
        save_parameters('SeisflowsParameters.json')
        save_paths('SeisflowsPaths.json')

        # submit job
        args = ('sbatch '
          + '--job-name=%s ' %  PAR.TITLE
          + '--output=%s ' % (PATH.SUBMIT+'/'+'output.log')
          + '--cpus-per-task=%d ' % PAR.NPROC
          + '--ntasks=%d ' % PAR.NTASK
          + '--time=%d ' % PAR.WALLTIME
          + findpath('system') +'/'+ 'slurm/wrapper_sbatch '
          + PATH.OUTPUT)

        subprocess.call(args, shell=1)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts
        """
        if PAR.VERBOSE >= 2:
            print 'running',funcname

        # save current state
        save_objects(join(PATH.OUTPUT,'SeisflowsObjects'))

        # save keyword arguments
        kwargspath = join(PATH.OUTPUT,'SeisflowsObjects',classname+'_kwargs')
        kwargsfile = join(kwargspath,funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile,kwargs)

        if hosts == 'all':
            # run on all available nodes
            args = ('srun '
              + '--wait=0 '
              + join(findpath('system'),'slurm/wrapper_srun ')
              + PATH.OUTPUT + ' '
              + classname + ' '
              + funcname)

        elif hosts == 'head':
            # run on head node
            args = ('srun '
              + '--wait=0 '
              + join(findpath('system'),'slurm/wrapper_srun_head ')
              + PATH.OUTPUT + ' '
              + classname + ' '
              + funcname)
        else:
            raise Exception

        subprocess.call(args, shell=1)


    def getnode(self):
        """ Gets number of running task
        """
        gid = os.getenv('SLURM_GTIDS').split(',')
        lid = int(os.getenv('SLURM_LOCALID'))
        return int(gid[lid])


    def mpiargs(self):
        return 'mpirun -np %d ' % PAR.NPROC
