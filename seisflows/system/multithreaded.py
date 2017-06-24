
import os
import sys
import numpy as np

from os.path import abspath, basename, join
from subprocess import Popen
from time import sleep

from seisflows.tools import unix
from seisflows.tools.tools import call, findpath, nproc, saveobj
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class multithreaded(custom_import('system', 'serial')):
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
        super(multithreaded, self).check()

        # number of tasks
        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        # number of cores per task
        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        # number of available cores
        if 'NPROCMAX' not in PAR:
            setattr(PAR, 'NPROCMAX', nproc())

        # maximum number of concurrent tasks
        if 'NTASKMAX' not in PAR:
            setattr(PAR, 'NTASKMAX', PAR.NPROCMAX/PAR.NPROC)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """ Runs tasks in serial or parallel on specified hosts
        """
        self.checkpoint()
        self.save_kwargs(classname, funcname, kwargs)

        if hosts == 'all':
            running_tasks = dict()
            queued_tasks = range(PAR.NTASK)

            # implements "work queue" pattern
            while queued_tasks or running_tasks:

                # launch queued tasks
                while len(queued_tasks) > 0 and \
                      len(running_tasks) < PAR.NTASKMAX:
                    i = queued_tasks.pop(0)
                    p = self._launch(classname, funcname, taskid=i)
                    running_tasks[i] = p
                    sleep(0.1)

                # checks status of running tasks
                for i, p in running_tasks.items():
                    if p.poll() != None:
                        running_tasks.pop(i)

                if running_tasks:
                    sleep(1.)

            print ''

        elif hosts == 'head':
            os.environ['SEISFLOWS_TASKID'] = str(0)
            func = getattr(__import__('seisflows_'+classname), funcname)
            func(**kwargs)

        else:
            raise KeyError('Bad keyword argument: system.run: hosts')


    ### private methods

    def _launch(self, classname, funcname, taskid=0):
        env = os.environ.copy().items()
        env += [['SEISFLOWS_TASKID', str(taskid)]]
        self.progress(taskid)

        p = Popen(
            findpath('seisflows.system') +'/'+ 'wrappers/run '
            + PATH.OUTPUT + ' '
            + classname + ' '
            + funcname,
            shell=True,
            env=dict(env))

        return p


    def save_kwargs(self, classname, funcname, kwargs):
        kwargspath = join(PATH.OUTPUT, 'kwargs')
        kwargsfile = join(kwargspath, classname+'_'+funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)

