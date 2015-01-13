
import os
import pickle
import sys
import subprocess
from os.path import abspath, join

from seisflows.tools import unix
from seisflows.tools.code import saveobj
from seisflows.tools.config import getmodule, findpath, ParameterObj

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

        if 'SYSTEM' not in PATH:
            setattr(PATH,'SYSTEM',join(PATH.GLOBAL,'system'))

        if 'SUBMIT' not in PATH:
            setattr(PATH,'SUBMIT',unix.pwd())

        if 'OUTPUT' not in PATH:
            setattr(PATH,'OUTPUT',join(PATH.SUBMIT,'output'))


    def submit(self, workflow):
        """Submits job
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)

        # save current state
        save_objects('SeisflowsObjects')
        save_parameters('SeisflowsParameters.json')
        save_paths('SeisflowsPaths.json')

        nodes = PAR.NTASK/16
        cores = PAR.NTASK%16
        hours = PAR.WALLTIME/60
        minutes = PAR.WALLTIME%60

        # construct resource list
        resources = 'walltime=%02d:%02d:00 ' % (hours, minutes)
        if nodes == 0:
            resources += ',nodes=1:ppn=%d' % cores
        elif cores == 0:
            resources += ',nodes=%d:ppn=16'% nodes
        else:
            resources += ',nodes=%d:ppn=16+1:ppn=%d' % (nodes, cores)

        args = ('qsub '
          + '-N %s ' % PAR.TITLE
          + '-o %s ' % (PATH.SUBMIT+'/'+'output.log')
          + '-l %s ' % resources
          + '-j %s ' % 'oe'
          + findpath('system') +'/'+ 'pbs/wrapper_qsub '
          + PATH.OUTPUT)
        print args # DEBUG

        subprocess.call(args, shell=1)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts
        """
        name = task.__name__

        if PAR.VERBOSE >= 2:
            print 'running',name

        # save current state
        save_objects(join(PATH.OUTPUT,'SeisflowsObjects'))

        # save keyword arguments
        kwargspath = join(PATH.OUTPUT,'SeisflowsObjects',classname+'_kwargs')
        kwargsfile = join(kwargspath,funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile,kwargs)

        if hosts == 'all':
            # run on all available nodes
            args = ('pbsdsh '
              + findpath('system') +'/'+ 'pbs/wrapper_pbsdsh '
              + PATH.OUTPUT + ' '
              + classname + ' '
              + funcname)
        elif hosts == 'head':
            # run on head node
            args = ('pbsdsh '
              + findpath('system') +'/'+ 'slurm/wrapper_pbsdsh '
              + PATH.OUTPUT + ' '
              + classname + ' '
              + funcname)
        else:
            raise Exception

        subprocess.call(args, shell=1)


    def getnode(self):
        """ Gets number of running task
        """
        return int(os.getenv('PBS_VNODENUM'))


    def mpiargs(self):
        return 'mpirun -np %d ' % PAR.NPROC
