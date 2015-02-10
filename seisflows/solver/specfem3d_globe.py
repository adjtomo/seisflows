
import subprocess
from glob import glob
from os.path import join

import numpy as np

import seisflows.seistools.specfem3d_globe as solvertools
from seisflows.seistools.io import load
from seisflows.seistools.shared import getpar, setpar

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class specfem3d_globe(loadclass('solver', 'base')):
    """ Python interface for SPECFEM3D_GLOBE

      See base class for method descriptions
    """

    if 0:
        # use isotropic model
        model_parameters = []
        model_parameters += ['reg1_rho']
        model_parameters += ['reg1_vp']
        model_parameters += ['reg1_vs']

        inversion_parameters = []
        inversion_parameters += ['reg1_rho']
        inversion_parameters += ['reg1_vp']
        inversion_parameters += ['reg1_vs']

        kernel_map = {
            'reg1_rho': 'reg1_rho_kernel',
            'reg1_vp': 'reg1_alpha_kernel',
            'reg1_vs': 'reg1_beta_kernel'}

    else:
        # use transversely isotropic model
        model_parameters = []
        model_parameters += ['reg1_rho']
        model_parameters += ['reg1_vpv']
        model_parameters += ['reg1_vph']
        model_parameters += ['reg1_vsv']
        model_parameters += ['reg1_vsh']
        model_parameters += ['reg1_eta']

        inversion_parameters = []
        inversion_parameters += ['reg1_rho']
        inversion_parameters += ['reg1_vpv']
        inversion_parameters += ['reg1_vph']
        inversion_parameters += ['reg1_vsv']
        inversion_parameters += ['reg1_vsh']
        inversion_parameters += ['reg1_eta']

        kernel_map = {
            'reg1_rho': 'reg1_rho_kernel',
            'reg1_eta': 'reg1_eta_kernel',
            'reg1_vph': 'reg1_alphah_kernel',
            'reg1_vpv': 'reg1_alphav_kernel',
            'reg1_vsv': 'reg1_betav_kernel',
            'reg1_vsh': 'reg1_betah_kernel'}


    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(specfem3d_globe, self).check()


    def generate_data(self, **model_kwargs):
        """ Generates data
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.getpath)
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xspecfem3D')

        unix.mv(self.data_path, 'traces/obs')
        self.export_traces(PATH.OUTPUT, 'traces/obs')


    def generate_mesh(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(model_type)

        self.initialize_solver_directories()
        unix.cd(self.getpath)

        if model_type == 'gll':
            assert (exists(model_path))
            unix.cp(glob(model_path +'/'+ '*'), self.model_databases)
            self.mpirun('bin/xmeshfem3D')
            self.export_model(PATH.OUTPUT +'/'+ model_name)

        else:
            raise NotImplementedError


    ### low-level solver interface

    def forward(self):
        """ Calls SPECFEM3D_GLOBE forward solver
        """
        solvertools.setpar('SIMULATION_TYPE', '1')
        solvertools.setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xspecfem3D')
        unix.mv(self.data_wildcard, 'traces/syn')


    def adjoint(self):
        """ Calls SPECFEM3D_GLOBE adjoint solver
        """
        solvertools.setpar('SIMULATION_TYPE', '3')
        solvertools.setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')
        self.mpirun('bin/xspecfem3D')


    ### postprocessing utilities

    def combine(self, path=''):
        """ combines SPECFEM3D_GLOBE kernels
        """
        dirs = unix.ls(path)

        # initialize kernels
        unix.cd(path)
        for key in self.model_parameters:
            if key not in self.inversion_parameters:
                for i in range(PAR.NPROC):
                    proc = '%06d' % i
                    name = self.kernel_map[key]
                    src = PATH.GLOBAL +'/'+ 'mesh' +'/'+ key +'/'+ proc
                    dst = path +'/'+ 'sum' +'/'+ 'proc'+proc+'_'+name+'.bin'
                    savebin(np.load(src), dst)

        # create temporary files and directories
        unix.cd(self.getpath)
        with open('kernels_list.txt', 'w') as file:
            file.write('\n'.join(dirs) + '\n')
        unix.mkdir('INPUT_KERNELS')
        unix.mkdir('OUTPUT_SUM')
        for dir in dirs:
            src = path +'/'+ dir
            dst = 'INPUT_KERNELS' +'/'+ dir
            unix.ln(src, dst)

        # sum kernels
        self.mpirun(PATH.SOLVER_BINARIES +'/'+ 'xsum_kernels')
        unix.mv('OUTPUT_SUM', path +'/'+ 'sum')

        # remove temporary files and directories
        unix.rm('INPUT_KERNELS')
        unix.rm('kernels_list.txt')

        unix.cd(path)


    def smooth(self, path='', tag='gradient', span=0.):
        """ smooths SPECFEM3D_GLOBE kernels
        """
        unix.cd(self.getpath)

        # list kernels
        kernels = []
        for name in self.model_parameters:
            if name in self.inversion_parameters:
                flag = True
            else:
                flag = False
            region, name = name.split('_')
            kernels = kernels + [[name, flag]]

        # smooth kernels
        for name, flag in kernels:
            if flag:
                print ' smoothing', name
                self.mpirun(
                    PATH.SOLVER_BINARIES +'/'+ 'xsmooth_sem '
                    + str(span) + ' '
                    + str(span) + ' '
                    + name + ' '
                    + path +'/'+ tag + '/ '
                    + self.model_databases + '/ ')

        # move kernels
        src = path +'/'+ tag
        dst = path +'/'+ '_nosmooth'
        unix.mkdir(dst)
        for name, flag in kernels:
            if flag:
                unix.mv(glob(src+'/*'+name+'.bin'), dst)
            else:
                unix.cp(glob(src+'/*'+name+'.bin'), dst)
        unix.rename('_smooth', '', glob(src+'/*'))
        print ''

        unix.cd(path)


    ### utility functions

    @property
    def data_wildcard(self):
        return glob('OUTPUT_FILES/*.sem.ascii')

    @property
    def model_databases(self):
        return join(self.getpath, 'OUTPUT_FILES/DATABASES_MPI')

    @property
    def source_prefix(self):
        return 'CMTSOLUTION'


