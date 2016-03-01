
import uuid

from os.path import abspath, join
from seisflows.tools import unix
from seisflows.tools.config import loadclass
from seisflows.tools.config import ParameterError, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class tiger_lg(loadclass('system', 'slurm_lg')):
    """ Specially designed system interface for tiger.princeton.edu

      See parent class for more information.
    """

    def check(self):
        """ Checks parameters and paths
        """

        if 'UUID' not in PAR:
            setattr(PAR, 'UUID', str(uuid.uuid4()))

        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH', join('/scratch/gpfs', unix.whoami(), 'seisflows', PAR.UUID))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', '')

        if 'NODESIZE' not in PAR:
            setattr(PAR, 'NODESIZE', 16)

        super(tiger_lg, self).check()


    def submit(self, *args, **kwargs):
        """Submits job
        """
        unix.ln(PATH.SCRATCH, PATH.SUBMIT + '/' + 'scratch')
        super(tiger_lg, self).submit(*args, **kwargs)

