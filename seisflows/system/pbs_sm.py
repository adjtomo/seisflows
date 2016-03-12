
from os.path import abspath, basename, join, dirname

from seisflows.tools import unix
from seisflows.tools.code import call, findpath, saveobj
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

      For important additional information, please see
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-configuration
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(pbs_sm, self).check()

        # check parameters
        if 'WALLTIME' not in PAR:
            setattr(PAR, 'WALLTIME', 30.)

        if 'MEMORY' not in PAR:
            setattr(PAR, 'MEMORY', 0)

        if 'NODESIZE' not in PAR:
            raise ParameterError(PAR, 'NODESIZE')

        if 'PBSARGS' not in PAR:
            setattr(PAR, 'PBSARGS', '')


    def submit(self, workflow):
        """Submits job
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)

        # save current state
        self.checkpoint()

        # construct resource list
        resources = []
        nodes = int(PAR.NTASK / PAR.NODESIZE)
        cores = PAR.NTASK % PAR.NODESIZE
        hours = int(PAR.WALLTIME / 60)
        minutes = PAR.WALLTIME % 60
        
        if PAR.WALLTIME:
            resources += ['walltime=%02d:%02d:00'%(hours, minutes)]
        if PAR.MEMORY:
            resources += ['mem=%dgb' % PAR.MEMORY]
        if nodes == 0:
            resources += ['nodes=1:ppn=%d'%(cores)]
        elif cores == 0:
            resources += ['nodes=%d:ppn=%d'%(nodes, PAR.NODESIZE)]
        else:
            resources += ['nodes=%d:ppn=%d+1:ppn=%d'%(nodes, PAR.NODESIZE, cores)]

        # construct arguments list
        call('qsub '
                + '%s ' % PAR.PBSARGS 
                + '-N %s '%PAR.TITLE
                + '-o %s '%(PATH.SUBMIT +'/'+ 'output.log')
                + '-l %s '%resources.join(',')
                + '-j %s '%'oe'
                + findpath('seisflows.system') +'/'+ 'wrappers/submit '
                + '-F %s '%PATH.OUTPUT)

