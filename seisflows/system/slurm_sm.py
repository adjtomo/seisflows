
import os
import subprocess
import sys
from os.path import abspath, join

from seisflows.tools import unix
from seisflows.tools.code import findpath, saveobj
from seisflows.tools.config import ParameterError, loadclass, \
    SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class slurm_sm(loadclass('system', 'mpi')):
    """ An interface through which to submit workflows, run tasks in serial or 
      parallel, and perform other system functions.

      By hiding environment details behind a python interface layer, these 
      classes provide a consistent command set across different computing
      environments.

      Intermediate files are written to a global scratch path PATH.SCRATCH,
      which must be accessible to all compute nodes.

      Optionally, users can provide a local scratch path PATH.LOCAL if each
      compute node has its own local filesystem.

      For more informations, see 
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-interfaces
    """


    def check(self):
        """ Checks parameters and paths
        """

        # check parameters
        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', unix.basename(abspath('.')))

        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        if 'SLURM_ARGS' not in PAR:
            setattr(PAR, 'SLURM_ARGS', '')

        # check paths
        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH', join(abspath('.'), 'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SUBMIT' not in PATH:
            setattr(PATH, 'SUBMIT', unix.pwd())

        if 'OUTPUT' not in PATH:
            setattr(PATH, 'OUTPUT', join(PATH.SUBMIT, 'output'))


    def submit(self, workflow):
        """ Submits workflow
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)

        self.checkpoint()

        # submit workflow
        unix.run('sbatch '
                + PAR.SLURM_ARGS + ' '
                + '--job-name=%s '%PAR.TITLE
                + '--output=%s '%(PATH.SUBMIT +'/'+ 'output.log')
                + '--cpus-per-task=%d '%PAR.NPROC
                + '--ntasks=%d '%PAR.NTASK
                + '--time=%d '%PAR.WALLTIME
                + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


