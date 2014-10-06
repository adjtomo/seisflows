
import subprocess
import sys

import numpy as np

from seisflows.tools import unix
from seisflows.tools.codetools import abspath, join
from seisflows.tools.configure import ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')


class serial(object):
    """ System interface class

      Provides an interface through which to submit jobs, run tasks in serial
      or parallel, and perform other system functions.

      One of several system interface classes included in SEISFLOWS in order
      to provide a consistent interface across different computer systems.
      Each class hides the details of submitting and running and jobs on one
      particular system.
    """


    def __init__(self):
        """ Class constructor
        """

        # check user supplied parameters
        if 'NTASK' not in PAR:
            setattr(PAR,'NTASK',1)

        if 'NPROC' not in PAR:
            setattr(PAR,'NPROC',1)

        if 'VERBOSE' not in PAR:
            setattr(PAR,'VERBOSE',1)

        if 'TITLE' not in PAR:
            setattr(PAR,'TITLE','')

        if 'SUBTITLE' not in PAR:
            setattr(PAR,'SUBTITLE','')

        # check user supplied paths
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
        """ Submits job
        """
        unix.mkdir(PATH.SYSTEM)

        workflow().main()


    def run(self,task,hosts='all',**kwargs):
        """ Runs tasks in serial or parallel on specified hosts
        """
        if hosts=='all':
            for itask in range(PAR.NTASK):
                self.setnode(itask)
                self.progress(itask)
                task(**kwargs)
            print ''

        elif hosts=='head':
            self.setnode(0)
            task(**kwargs)

        else:
            task(**kwargs)


    def getnode(self):
        "Gets number of running task"
        return int(np.loadtxt(PATH.SYSTEM+'/'+'nodenum'))


    def setnode(self,itask):
        "Sets number of running task"
        np.savetxt(PATH.SYSTEM+'/'+'nodenum',[itask])


    def mpiexec(self):
        "Wrapper for mpiexec"
        return 'mpiexec -np %d ' % PAR.NPROC


    def progress(self,itask=None):
        "Prints status updates"
        if PAR.VERBOSE and PAR.NTASK > 1:
            print ' task '+'%02d'%(itask+1)+' of '+'%02d'%PAR.NTASK
