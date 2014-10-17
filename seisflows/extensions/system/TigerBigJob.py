
from seisflows.tools import unix
from seisflows.tools.codetools import abspath, join
from seisflows.tools.configtools import getclass, GlobalStruct

PAR = GlobalStruct('parameters')
PATH = GlobalStruct('paths')


class TigerBigJob(getclass('system','slurm_big_job')):

    def __init__(self):
        """ Class constructor
        """

        # check parameters
        if 'NTASK' not in PAR:
            raise Exception

        if 'NPROC' not in PAR:
            raise Exception

        if 'CPUS_PER_NODE' not in PAR:
            setattr(PAR,'CPUS_PER_NODE',16)

        if 'WALLTIME' not in PAR:
            setattr(PAR,'WALLTIME',30.)

        if 'STEPTIME' not in PAR:
            setattr(PAR,'STEPTIME',30.)

        if 'VERBOSE' not in PAR:
            setattr(PAR,'VERBOSE',1)

        if 'TITLE' not in PAR:
            setattr(PAR,'TITLE',unix.basename(abspath('.')))

        if 'SUBTITLE' not in PAR:
            setattr(PAR,'SUBTITLE',unix.basename(abspath('..')))

        if 'RETRY' not in PAR:
            PAR.RETRY = False

        if 'SLEEP' not in PAR:
            PAR.SLEEP = 30.

        # check paths
        if 'GLOBAL' not in PATH:
            setattr(PATH,'GLOBAL',join('/scratch/gpfs',unix.whoami(), \
                PAR.SUBTITLE,PAR.TITLE))

        if 'LOCAL' not in PATH:
            setattr(PATH,'LOCAL','')

        if 'SYSTEM' not in PATH:
            setattr(PATH,'SYSTEM',join(PATH.GLOBAL,'system'))

        if 'SUBMIT' not in PATH:
            setattr(PATH,'SUBMIT',unix.pwd())


    def submit(self,*args,**kwargs):
        """Submits job
        """
        unix.ln(PATH.GLOBAL,PATH.SUBMIT+'/'+'scratch')
        super(TigerBigJob,self).submit(*args,**kwargs)
