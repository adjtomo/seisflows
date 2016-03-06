
from getpass import getuser
from os.path import abspath, join, exists
from uuid import uuid4

from seisflows.tools import unix
from seisflows.tools.config import ParameterError, SeisflowsParameters, SeisflowsPaths, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class tiger_md(custom_import('system', 'slurm_md')):
    """ Specially designed system interface for tiger.princeton.edu

      See parent class for more information.
    """

    def check(self):
        """ Checks parameters and paths
        """

        if 'UUID' not in PAR:
            setattr(PAR, 'UUID', str(uuid4()))
 
        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH', join('/scratch/gpfs', getuser(), 'seisflows', PAR.UUID))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', '')

        super(tiger_md, self).check()


    def submit(self, *args, **kwargs):
        """ Submits job
        """
        if not exists(PATH.SUBMIT + '/' + 'scratch'):
            unix.ln(PATH.SCRATCH, PATH.SUBMIT + '/' + 'scratch')

        super(tiger_md, self).submit(*args, **kwargs)
