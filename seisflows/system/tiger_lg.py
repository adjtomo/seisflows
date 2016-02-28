
from os.path import abspath, join
from seisflows.tools import unix
from seisflows.tools.config import loadclass
from seisflows.tools.config import ParameterError, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class tiger_lg(loadclass('system', 'slurm_lg')):
    """ Specially designed system interface for tiger.princeton.edu

      For more informations, see 
      http://seisflows.readthedocs.org/en/latest/manual/manual.html#system-interfaces
    """

    def check(self):
        """ Checks parameters and paths
        """

        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', unix.basename(abspath('..')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR, 'SUBTITLE', unix.basename(abspath('.')))

        if 'SCRATCH' not in PATH:
            setattr(PATH, 'SCRATCH',
                    join('/scratch/gpfs', unix.whoami(), PAR.TITLE, PAR.SUBTITLE))

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

