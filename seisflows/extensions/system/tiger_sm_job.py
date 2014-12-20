from seisflows.tools import unix
from seisflows.tools.code import abspath, exists, join
from seisflows.tools.config import loadclass, ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

save_objects = OBJ.save
save_parameters = PAR.save
save_paths = PATH.save


class tiger_sm_job(loadclass('system', 'slurm_sm_job')):
    def check(self):
        """ Checks parameters and paths
        """

        if 'TITLE' not in PAR:
            setattr(PAR, 'TITLE', unix.basename(abspath('.')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR, 'SUBTITLE', unix.basename(abspath('..')))

        if 'SUBDIRS' not in PATH:
            setattr(PATH, 'SUBDIRS', join(PAR.SUBTITLE, PAR.TITLE))

        if 'GLOBAL' not in PATH:
            setattr(PATH, 'GLOBAL',
                    join('/scratch/gpfs', unix.whoami(), PATH.SUBDIRS))

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', '')

        super(self.__class__, self).check()


    def submit(self, *args, **kwargs):
        """Submits job
        """
        if not exists(PATH.SUBMIT + '/' + 'scratch'):
            unix.ln(PATH.GLOBAL, PATH.SUBMIT + '/' + 'scratch')

        super(self.__class__, self).submit(*args, **kwargs)
