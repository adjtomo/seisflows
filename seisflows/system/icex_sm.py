
from os.path import abspath, join
from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass
from seisflows.tools.config import ParameterError, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class icex_sm(loadclass('system', 'lsf_sm')):
    """ Specially designed system interface for ICEXDEV

      By hiding environment details behind a python interface layer, these 
      classes provide a consistent command set across different computing
      environments.

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

        if 'GLOBAL' not in PATH:
            setattr(PATH, 'GLOBAL',
                    join('/scratch/gpfs', unix.whoami(), PAR.TITLE, PAR.SUBTITLE))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', '')

        super(icex_sm, self).check()


    def submit(self, *args, **kwargs):
        """Submits job
        """
        if not exists(PATH.SUBMIT + '/' + 'scratch'):
            unix.ln(PATH.GLOBAL, PATH.SUBMIT + '/' + 'scratch')

        super(icex_sm, self).submit(*args, **kwargs)
