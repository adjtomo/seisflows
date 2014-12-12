
import os
import pickle
import sys
import subprocess

from seisflows.tools import unix
from seisflows.tools.code import abspath, join, saveobj
from seisflows.tools.config import getmodule, getpath, ParameterObj

PAR = ParameterObj('parameters')
PATH = ParameterObj('paths')


class pbs_torque(object):
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

        # store parameters
        saveobj(join(PATH.SUBMIT_DIR,'parameters.p'),PAR.vars)
        saveobj(join(PATH.SUBMIT_DIR,'paths.p'),PATH.vars)

        # construct resource list
        nodes = PAR.NTASK/16
        cores = PAR.NTASK%16
        if nodes == 0:
            rlist = 'nodes=1:ppn=%d' % cores
        elif cores == 0:
            rlist = 'nodes=%d:ppn=16'% nodes
        else:
            rlist = 'nodes=%d:ppn=16+1:ppn=%d' % (nodes,cores)

        hh = PAR.WALLTIME/60
        mm = PAR.WALLTIME%60
        resources = "-l %s,walltime=%02d:%02d:00 " % (rlist,hh,mm)

        # construct variables list
        variables = ''

        args = ('qsub '
          + resources
          + variables
          + '-N %s ' %  PAR.TITLE
          + '-o %s ' % (PATH.SUBMIT_DIR+'/'+'output.log')
          + '-j %s ' % 'oe'
          + getpath('system') +'/'+ 'pbs/wrapper_qsub '
          + PATH.SUBMIT_DIR + ' '
          + getmodule(workflow))

        print args

        subprocess.call(args, shell=True)


    def run(self,task,hosts='all',**kwargs):
        """  Runs tasks in serial or parallel on specified hosts
        """
        name = task.__name__

        if PAR.VERBOSE >= 2:
            print 'running',name

        # store function arguments
        unix.mkdir(PATH.SYSTEM)
        file = PATH.SYSTEM+'/'+name+'.p'
        saveobj(file,kwargs)

        if hosts == 'all':
            # run on all available nodes
            args = ('pbsdsh '
              + getpath('system') +'/'+ 'pbs/wrapper_pbsdsh '
              + PATH.SUBMIT_DIR + ' '
              + getmodule(task) + ' '
              + name)
        elif hosts == 'head':
            # run on head node
            args = ('pbsdsh '
              + getpath('system') +'/'+ 'slurm/wrapper_pbsdsh '
              + PATH.SUBMIT_DIR + ' '
              + getmodule(task) + ' '
              + name)
        else:
            raise Exception

        subprocess.call(args, shell=True)


    def getnode(self):
        """ Gets number of running task
        """
        return int(os.getenv('PBS_VNODENUM'))


    def mpiargs(self):
        return 'mpirun -np %d ' % PAR.NPROC
