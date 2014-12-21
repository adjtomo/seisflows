import subprocess

from seisflows.tools import unix
from seisflows.tools.code import glob
from seisflows.tools.config import loadclass, ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import preprocess


class specfem3d_legacy(loadclass('solver', 'specfem3d')):

    def mpirun(self, runfile, args='', outfile='/dev/null'):
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

    def combine(self, path=''):
        """ combines SPECFEM3D kernels
        """
        unix.cd(self.spath)

        # create temporary files and directories needed by xsum_kernels
        dirs = unix.ls(path)
        with open('../kernels_list.txt', 'w') as file:
            file.write('\n'.join(dirs) + '\n')
        unix.mkdir('INPUT_KERNELS')
        unix.mkdir('OUTPUT_SUM')
        for dir in dirs:
            src = path + '/' + dir
            dst = unix.pwd() + '/' + 'INPUT_KERNELS' + '/' + dir
            unix.ln(src, dst)

        # sum kernels
        self.mpirun(PATH.SOLVER_BINARIES + '/' + 'xsum_kernels')
        unix.mv('OUTPUT_SUM', path + '/' + 'sum')

        # remove temporary files and directories
        unix.rm('INPUT_KERNELS')
        unix.rm('../kernels_list.txt')

        unix.cd(path)

    def smooth(self, path='', tag='grad', span=0.):
        """ smooths SPECFEM3D kernels
        """
        unix.cd(self.spath)

        # list kernels
        kernels = []
        for name in self.model_parameters:
            if name in self.inversion_parameters:
                flag = True
            else:
                flag = False
            kernels = kernels + [[name, flag]]

        # smooth kernels
        for name, flag in kernels:
            if flag:
                print ' smoothing', name
                self.mpirun(
                    PATH.SOLVER_BINARIES + '/' + 'xsmooth_vol_data ',
                    name + ' '
                    + path + '/' + tag + '_nosmooth/' + ' '
                    + path + '/' + tag + '/' + ' '
                    + str(span) + ' '
                    + str(span) + ' ' + '1')

        # move kernels
        src = path +'/'+ tag
        dst = path +'/'+ tag + '_nosmooth'
        unix.mkdir(dst)
        for name, flag in kernels:
            if flag:
                unix.mv(glob(src+'/*'+name+'.bin'), dst)
            else:
                unix.cp(glob(src+'/*'+name+'.bin'), dst)
        unix.rename('_smooth', '', glob(src+'/*'))
        print ''


        unix.cd(path)

