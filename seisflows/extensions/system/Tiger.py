
from seisflows.tools import unix
from seisflows.tools.codetools import abspath, exists, join
from seisflows.tools.configtools import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class Tiger(loadclass('system','slurm')):

    def check(self):
        """ Class constructor
        """

        if 'TITLE' not in PAR:
            setattr(PAR,'TITLE',unix.basename(abspath('.')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR,'SUBTITLE',unix.basename(abspath('..')))

        if 'UID' not in PATH:
            setattr(PATH,'UID',join(PAR.SUBTITLE,PAR.TITLE))

        # check parameters
        if 'NTASK' not in PAR:
            raise Exception

        if 'NPROC' not in PAR:
            raise Exception

        if 'WALLTIME' not in PAR:
            setattr(PAR,'WALLTIME',30.)

        if 'VERBOSE' not in PAR:
            setattr(PAR,'VERBOSE',1)

        # check paths
        if 'GLOBAL' not in PATH:
            setattr(PATH,'GLOBAL',join('/scratch/gpfs',unix.whoami(),PATH.UID))

        if 'LOCAL' not in PATH:
            setattr(PATH,'LOCAL','')

        if 'SYSTEM' not in PATH:
            setattr(PATH,'SYSTEM',join(PATH.GLOBAL,'system'))

        if 'SUBMIT' not in PATH:
            setattr(PATH,'SUBMIT',unix.pwd())


    def submit(self,*args,**kwargs):
        """Submits job
        """
        if not exists(PATH.SUBMIT+'/'+'scratch'):
            unix.ln(PATH.GLOBAL,PATH.SUBMIT+'/'+'scratch')

        super(Tiger,self).submit(*args,**kwargs)
