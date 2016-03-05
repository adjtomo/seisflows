import os
import subprocess
from os.path import abspath, basename, join, dirname

from seisflows.tools import unix
from seisflows.tools.code import findpath, saveobj
from seisflows.tools.config import ParameterError, custom_import, \
    SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class pbs_sm(custom_import('system', 'mpi')):
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
            setattr(PAR, 'TITLE', basename(abspath('.')))

        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        if 'MEMORY' not in PAR:
            raise ParameterError(PAR, 'MEMORY')

        if 'VERBOSE' not in PAR:
            setattr(PAR, 'VERBOSE', 1)

        if 'NTASK' not in PAR:
            raise ParameterError(PAR, 'NTASK')

        if 'NPROC' not in PAR:
            raise ParameterError(PAR, 'NPROC')

        if 'NODESIZE' not in PAR:
            raise ParameterError(PAR, 'NODESIZE')

        # check paths
        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH', join(abspath('.'), 'scratch'))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'SYSTEM' not in PATH:
            setattr(PATH, 'SYSTEM', join(PATH.SCRATCH, 'system'))

        if 'SUBMIT' not in PATH:
            setattr(PATH, 'SUBMIT', abspath('.'))

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
        nodes = int(PAR.NTASK / PAR.NODESIZE)
        cores = PAR.NTASK % PAR.NODESIZE
        hours = int(PAR.WALLTIME / 60)
        minutes = PAR.WALLTIME % 60
        resources = 'walltime=%02d:%02d:00'%(hours, minutes)
        if nodes == 0:
            resources += ',mem=%dgb,nodes=1:ppn=%d'%(PAR.MEMORY, cores)
        elif cores == 0:
            resources += ',mem=%dgb,nodes=%d:ppn=%d'%(PAR.MEMORY, nodes, PAR.NODESIZE)
        else:
            resources += ',mem=%dgb,nodes=%d:ppn=%d+1:ppn=%d'%(PAR.MEMORY, nodes, PAR.NODESIZE, cores)

        # construct arguments list
        unix.run('qsub '
                + '-N %s '%PAR.TITLE
                + '-o %s '%(PATH.SUBMIT +'/'+ 'output.log')
                + '-l %s '%resources
                + '-j %s '%'oe'
                + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + '-F %s '%PATH.OUTPUT)

