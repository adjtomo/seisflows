
from os.path import abspath, join
from seisflows.tools import unix
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class tiger_lg(loadclass('system', 'slurm_lg')):
    """ Specially designed system interface for tiger.princeton.edu

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

        if 'NPROC_PER_NODE' not in PAR:
            setattr(PAR, 'NPROC_PER_NODE', 16)

        super(tiger_lg, self).check()


    def submit(self, *args, **kwargs):
        """Submits job
        """
        unix.ln(PATH.GLOBAL, PATH.SUBMIT + '/' + 'scratch')
        super(tiger_lg, self).submit(*args, **kwargs)

