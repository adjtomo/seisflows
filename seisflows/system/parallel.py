from os.path import abspath, join

import numpy as np

from seisflows.tools import unix
from seisflows.tools.config import ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


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

        if 'NPROC_MAX' not in PAR:
            raise Exception


    def run(self, classname, funcname, hosts='all', **kwargs):
        """ Runs tasks in serial or parallel on specified hosts
        """
        unix.mkdir(PATH.SYSTEM)

        if hosts == 'all':
            itask = 0
            running_tasks = []
            queued_tasks = range(min(PAR.NPROC_MAX/PAR.NPROC, PAR.NTASK))

            while 1:
                # check running tasks
                for i, p in running_tasks:
                    if p.poll():
                        itask += 1
                        running_tasks.pop(i)
                        queued_tasks.append(itask)

                if itask == PAR.NTASK-1:
                    break

                # launch new tasks
                for i in queued_tasks:
                    p = self.launch(classname, funcname, itask=i)
                    running_tasks.add(i, p)
            print ''

        elif hosts == 'head':
            self.setnode(0)
            func = getattr(__import__(classname), funcname)
            func(**kwargs)

        else:
            task(**kwargs)


    def launch(classname, funcname, itask=0):
        self.progress(itask)

        env = os.environ.copy().items()
        env += ['SEIFLOWS_TASKID', itask]

        p = subprocess.call(myscript
            + classname + ' '
            + funcname + ' '
            shell=True,
            env=dict(env))

        return p



