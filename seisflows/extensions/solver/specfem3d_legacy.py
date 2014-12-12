
import subprocess

from seisflows.tools import unix
from seisflows.tools.code import glob
from seisflows.tools.config import loadclass, ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class specfem3d_legacy(loadclass('solver','specfem3d')):

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(specfem3d_legacy,self).check()

        if 'system' not in OBJ:
            raise Exception("Undefined Exception")

        global system
        import system


    def mpirun(self,runfile,args='',outfile='/dev/null'):
        """ Wrapper for mpirun
        """
        unix.cd('bin')

        with open(outfile) as f:
            subprocess.call(
                  system.mpiargs() +
                      unix.basename(runfile) +
                      args,
                  shell=True,
                  stdout=f)
        unix.cd('..')


    def combine(self,path=''):
        """ combines SPECFEM3D kernels
        """
        unix.cd(self.path)

        # create temporary files and directories needed by xsum_kernels
        dirs = unix.ls(path)
        with open('../kernels_list.txt','w') as file:
            file.write('\n'.join(dirs)+'\n')
        unix.mkdir('INPUT_KERNELS')
        unix.mkdir('OUTPUT_SUM')
        for dir in dirs:
            src = path + '/' + dir
            dst = unix.pwd() +'/'+ 'INPUT_KERNELS' +'/'+ dir
            unix.ln(src,dst)

        # sum kernels
        self.mpirun(PATH.SOLVER_BINARIES+'/'+'xsum_kernels')
        unix.mv('OUTPUT_SUM',path+'/'+'sum')

        # remove temporary files and directories
        unix.rm('INPUT_KERNELS')
        unix.rm('../kernels_list.txt')

        unix.cd(path)


    def smooth(self,path='',tag='grad',span=0):
        """ smooths SPECFEM3D kernels
        """
        unix.cd(self.path)

        unix.mv(path+'/'+tag,path+'/'+tag+'_nosmooth')
        unix.mkdir(path+'/'+tag)

        # prepare list
        kernel_list = []
        for key in self.model_parameters:
            if key in self.inversion_parameters:
                smoothing_flag = True
            else:
                smoothing_flag = False
            kernel_list = kernel_list + [[key,smoothing_flag]]

        for kernel_name,smoothing_flag in kernel_list:
            if smoothing_flag:
                # run smoothing
                print ' smoothing', kernel_name
                self.mpirun(
                  PATH.SOLVER_BINARIES+'/'+'xsmooth_vol_data ',
                  kernel_name + ' '
                  + path + '/' + tag+'_nosmooth/' + ' '
                  + path + '/' + tag+'/' + ' '
                  + str(span) + ' '
                  + str(span) + ' ' + '1')
            else:
                src = glob(path+'/'+tag+'_nosmooth/*'+kernel_name+'.bin')
                dst = path+'/'+tag+'/'
                unix.cp(src,dst)

        unix.rename('_smooth','',glob(path+'/'+tag+'/*_smooth.bin'))
        print ''

        unix.cd(path)

