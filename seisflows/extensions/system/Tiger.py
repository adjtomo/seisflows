
from seisflows.tools import unix
from seisflows.tools.codetools import abspath, exists, join
from seisflows.tools.configtools import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')


class Tiger(getclass('system','slurm')):

    def __init__(self):
        """ Class constructor
        """

        # check parameters
        if 'NTASK' not in PAR:
            raise Exception

        if 'NPROC' not in PAR:
            raise Exception

        if 'WALLTIME' not in PAR:
            setattr(PAR,'WALLTIME',30.)

        if 'VERBOSE' not in PAR:
            setattr(PAR,'VERBOSE',1)

        if 'TITLE' not in PAR:
            setattr(PAR,'TITLE',unix.basename(abspath('.')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR,'SUBTITLE',unix.basename(abspath('..')))

        # check paths
        if 'APPEND_TO_PATH' not in PATH:
            setattr(PATH,'APPEND_TO_PATH',join(PAR.SUBTITLE,PAR.TITLE))

        if 'GLOBAL' not in PATH:
            setattr(PATH,'GLOBAL',join('/scratch/gpfs',unix.whoami(), \
                PATH.APPEND_TO_PATH))

        if 'LOCAL' not in PATH:
            setattr(PATH,'LOCAL','')

        if 'SYSTEM' not in PATH:
            setattr(PATH,'SYSTEM',join(PATH.GLOBAL,'system'))

        if 'SUBMIT_DIR' not in PATH:
            setattr(PATH,'SUBMIT_DIR',unix.pwd())

        if 'SUBMIT_HOST' not in PATH:
            setattr(PATH,'SUBMIT_HOST',unix.hostname())


    def submit(self,*args,**kwargs):
        """Submits job
        """
        if not exists(PATH.SUBMIT_DIR+'/'+'scratch'):
            unix.ln(PATH.GLOBAL,PATH.SUBMIT_DIR+'/'+'scratch')

        super(Tiger,self).submit(*args,**kwargs)
