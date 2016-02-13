import os
import subprocess
from os.path import abspath, join

from seisflows.tools import unix
from seisflows.tools.code import saveobj
from seisflows.tools.config import ParameterError, findpath, loadclass, \
    SeisflowsObjects, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class pbs_sm(object):
    """ An interface through which to submit workflows, run tasks in serial or 
      parallel, and perform other system functions.

      By hiding environment details behind a python interface layer, these 
      classes provide a consistent command set across different computing
      environments.

      For more informations, see 
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-interfaces
    """

    def __init__(self):
        """ Checks parameters and paths
        """

        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', unix.basename(abspath('..')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR, 'SUBTITLE', unix.basename(abspath('.')))

        # check parameters
        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        if 'NODESIZE' not in PAR:
            raise ParameterError(PAR, 'NODESIZE')

        # check paths
        if 'GLOBAL' not in PATH:
            setattr(PATH, 'GLOBAL', join(abspath('.'), 'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SYSTEM' not in PATH:
            setattr(PATH, 'SYSTEM', join(PATH.GLOBAL, 'system'))

        if 'SUBMIT' not in PATH:
            setattr(PATH, 'SUBMIT', unix.pwd())

        if 'OUTPUT' not in PATH:
            setattr(PATH, 'OUTPUT', join(PATH.SUBMIT, 'output'))


    def submit(self, workflow):
        """Submits job
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)

        # save current state
        self.checkpoint()

        # construct resource list
        nodes = PAR.NTASK/PAR.NODESIZE
        cores = PAR.NTASK%PAR.NODESIZE
        hours = PAR.WALLTIME/60
        minutes = PAR.WALLTIME%60
        resources = 'walltime=%02d:%02d:00 '%(hours, minutes)
        if nodes == 0:
            resources += ',nodes=1:ppn=%d'%cores
        elif cores == 0:
            resources += ',nodes=%d:ppn=16'%nodes
        else:
            resources += ',nodes=%d:ppn=16+1:ppn=%d'%(nodes, cores)

        # construct arguments list
        unix.run('qsub '
                + '-N %s '%PAR.TITLE
                + '-o %s '%(PATH.SUBMIT +'/'+ 'output.log')
                + '-l %s '%resources
                + '-j %s '%'oe'
                + findpath('system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT,
                shell=True)


    def run(self, classname, funcname, hosts='all', **kwargs):
        """  Runs tasks in serial or parallel on specified hosts
        """
        self.checkpoint()
        self.save_kwargs(classname, funcname, kwargs)

        if hosts == 'all':
            # run on all available nodes
            unix.run('srun '
                    + '--wait=0 '
                    + join(findpath('system'), 'wrappers/run_pbsdsh ')
                    + PATH.OUTPUT + ' '
                    + classname + ' '
                    + funcname)

        elif hosts == 'head':
            # run on head node
            unix.run('srun '
                    + '--wait=0 '
                    + join(findpath('system'), 'wrappers/run_pbsdsh_head ')
                    + PATH.OUTPUT + ' '
                    + classname + ' '
                    + funcname + ' '
                    + findpath('seisflows'),
                    shell=True)


    def getnode(self):
        """ Gets number of running task
        """
        return int(os.getenv('PBS_VNODENUM'))

    def mpiargs(self):
        return 'mpirun -np %d '%PAR.NPROC

    def save_kwargs(self, classname, funcname, kwargs):
        kwargspath = join(PATH.OUTPUT, 'SeisflowsObjects', classname+'_kwargs')
        kwargsfile = join(kwargspath, funcname+'.p')
        unix.mkdir(kwargspath)
        saveobj(kwargsfile, kwargs)

