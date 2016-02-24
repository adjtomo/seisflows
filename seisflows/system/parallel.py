
import os

from os.path import abspath, join
from subprocess import Popen
from time import sleep

import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import saveobj
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, findpath, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class parallel(loadclass('system', 'serial')):
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
        super(parallel, self).check()

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

            while True:
                # check running tasks
                for i, p in running_tasks.items():
                    if p.poll() != None:
                        running_tasks.pop(i)

                # launch new tasks
                while len(running_tasks) < PAR.NPROCMAX and queued_tasks:
                    i = queued_tasks.pop(0)
                    p = self._launch(classname, funcname, itask=i)
                    running_tasks[i] = p

                if running_tasks:
                    sleep(1)
                    continue

                if not queued_tasks:
                    break

            print ''

        elif hosts == 'head':
            self.setnode(0)
            func = getattr(__import__(classname), funcname)
            func(**kwargs)

        else:
            task(**kwargs)


    ### private methods

    def _launch(self, classname, funcname, itask=0):
        self.progress(itask)

        env = os.environ.copy().items()
        env += [['SEISFLOWS_TASKID', str(itask)]]

        p = Popen(
            findpath('system') +'/'+ 'wrappers/run '
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

