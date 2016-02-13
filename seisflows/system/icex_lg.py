
from os.path import abspath, join
from seisflows.tools import unix
from seisflows.tools.config import loadclass
from seisflows.tools.config import ParameterError, SeisflowsParameters, SeisflowsPaths

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class icex_lg(loadclass('system', 'lsf_lg')):
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

        if 'NODESIZE' not in PAR:
            setattr(PAR, 'NODESIZE', 16)

        super(icex_lg, self).check()


    def submit(self, workflow):
        """ Submits workflow
        """
        unix.mkdir(PATH.OUTPUT)
        unix.cd(PATH.OUTPUT)
        unix.mkdir(PATH.SUBMIT+'/'+'output.lsf')

        self.save_objects()
        self.save_parameters()
        self.save_paths()

        # prepare bsub arguments
        unix.run('bsub '
                + '-a intelmpi '
                + '-J %s ' % PAR.SUBTITLE
                + '-q LAURE_USERS '
                + '-o %s ' % (PATH.SUBMIT+'/'+'output.log')
                + '-n %d ' % 16
                + '-e %s ' % (PATH.SUBMIT+'/'+'error.log')
                + '-R "span[ptile=%d' % PAR.NODESIZE + ']" '
                + '-W %d:00 ' % PAR.WALLTIME
                +  findpath('system') +'/'+ 'wrappers/submit '
                + PATH.OUTPUT)


    def mpiargs(self):
        #return 'mpirun '
        #return ('/apps/lsf/cluster_ICEX/8.3/linux2.6-glibc2.3-x86_64/bin/mpirun.lsf '
        #        + '-genv I_MPI_EXTRA_FILESYSTEM 1 -genv I_MPI_EXTRA_FILESYSTEM_LIST lustre '
        #        + '-genv I_MPI_PIN 0 -genv I_MPI_FALLBACK 0 -_MSG_SIZE 4194304 -pam ')
        return ('/apps/lsf/cluster_ICEX/8.3/linux2.6-glibc2.3-x86_64/bin/mpirun.lsf '
                + '-genv I_MPI_EXTRA_FILESYSTEM 1 -genv I_MPI_EXTRA_FILESYSTEM_LIST lustre '
                + '-genv I_MPI_PIN 0 -genv I_MPI_FALLBACK 0 -genv I_MPI_RDMA_RNDV_WRITE 1 -genv I_MPI_RDMA_MAX_MSG_SIZE 4194304 -pam '
                + ' "-n %s " ' % PAR.NPROC )

