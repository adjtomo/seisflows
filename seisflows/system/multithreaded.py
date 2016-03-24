
import os

from os.path import abspath, basename, join
from subprocess import Popen
from time import sleep

import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import call, findpath, saveobj
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


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

        if 'NPROCMAX' not in PAR:
            raise Exception


    def run(self, classname, funcname, hosts='all', **kwargs):
        """ Runs tasks in serial or parallel on specified hosts
        """
        unix.mkdir(PATH.SYSTEM)

        self.checkpoint()
        self.save_kwargs(classname, funcname, kwargs)

        if hosts == 'all':
            running_tasks = dict()
            queued_tasks = range(PAR.NTASK)

            # implements "work queue" pattern
            while queued_tasks or running_tasks:

                # launch queued tasks
                while len(queued_tasks) > 0 and \
                      len(running_tasks) < PAR.NPROCMAX:
                    i = queued_tasks.pop(0)
                    p = self._launch(classname, funcname, itask=i)
                    running_tasks[i] = p

                # checks status of running tasks
                for i, p in running_tasks.items():
                    if p.poll() != None:
                        running_tasks.pop(i)

                if running_tasks:
                    sleep(1)

            print ''

        elif hosts == 'head':
            self.setnode(0)
            func = getattr(__import__(classname), funcname)
            func(**kwargs)

        else:
            raise(KeyError('Hosts parameter not set/recognized.'))


    ### private methods

    def _launch(self, classname, funcname, itask=0):
        self.progress(itask)

        env = os.environ.copy().items()
        env += [['SEISFLOWS_TASKID', str(itask)]]

        p = Popen(
            findpath('seisflows.system') +'/'+ 'wrappers/run '
            + PATH.OUTPUT + ' '
            + classname + ' '
            + funcname,
            shell=True,
            env=dict(env))

        return p


    def save_kwargs(self, classname, funcname, kwargs):
        kwargspath = join(PATH.OUTPUT, 'SeisflowsObjects', classname+'_kwargs')
        kwargsfile = join(kwargspath, funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)

