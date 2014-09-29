
import os
import pickle
import sys
import subprocess

from seisflows.tools import unix
from seisflows.tools.codetools import abspath, join, saveobj
from seisflows.tools.configure import getmodule, getpath, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')


class slurm(object):
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

      if 'SYSTEM' not in PATH:
          setattr(PATH,'SYSTEM',join(PATH.GLOBAL,'system'))

      if 'SUBMIT_DIR' not in PATH:
          setattr(PATH,'SUBMIT_DIR',unix.pwd())

      if 'SUBMIT_HOST' not in PATH:
          setattr(PATH,'SUBMIT_HOST',unix.hostname())


  def submit(self,workflow):
      """Submits job
      """
      unix.cd(PATH.SUBMIT_DIR)
      unix.rm(PATH.GLOBAL)
      unix.mkdir(PATH.GLOBAL)

      # store parameters
      unix.mkdir(PATH.SYSTEM)
      saveobj(join(PATH.SYSTEM,'parameters.p'),
          PAR.vars)
      saveobj(join(PATH.SYSTEM,'paths.p'),
          PATH.vars)

      args = ('sbatch '
        + '--job-name=%s ' %  PAR.TITLE
        + '--output=%s ' % (PATH.SUBMIT_DIR+'/'+'output.log')
        + '--cpus-per-task=%d ' % PAR.NPROC
        + '--ntasks=%d ' % PAR.NTASK
        + '--time=%d ' % PAR.WALLTIME
        + getpath('system') +'/'+ 'slurm/wrapper_sbatch '
        + PATH.SYSTEM + ' '
        + getmodule(workflow))

      subprocess.call(args, shell=1)
   

  def run(self,task,hosts='all',**kwargs):
      """  Runs tasks in serial or parallel on specified hosts
      """
      name = task.__name__

      if PAR.VERBOSE >= 2:
        print 'running',name

      # store function arguments
      file = PATH.SYSTEM+'/'+name+'.p'
      saveobj(file,kwargs)

      if hosts == 'all':
        # run on all available nodes
        args = ('srun '
          + '--wait=0 '
          + getpath('system') +'/'+ 'slurm/wrapper_srun '
          + PATH.SYSTEM + ' '
          + getmodule(task) + ' '
          + name)
      elif hosts == 'head':
        # run on head node
        args = ('srun '
          + '--wait=0 '
          + getpath('system') +'/'+ 'slurm/wrapper_srun_head '
          + PATH.SYSTEM + ' '
          + getmodule(task) + ' '
          + name)
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

