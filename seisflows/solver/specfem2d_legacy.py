
import subprocess
from os.path import join
from glob import glob

import numpy as np

import seisflows.seistools.specfem2d as solvertools
from seisflows.seistools.shared import getpar, setpar
from seisflows.seistools.io import splitvec

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import preprocess


class specfem2d_legacy(loadclass('solver', 'base')):
    """ Python interface for SPECFEM2D

      See base class for method descriptions
    """

    if PAR.MATERIALS == 'Acoustic':
        # SPECFEM2D consists of separate acoustic, elastic and anisotropic 
        # solvers. Prior to the addition of the acoustic solver, acoustic 
        # simulations could be carried out through the elastic solver as SH-
        # simulations, the mathematics being essentially the same.
        parameters = []
        parameters += ['vs']


    def check(self):
        """ Checks parameters and paths
        """
        super(specfem2d_legacy, self).check()

        # check time stepping parameters
        if 'NT' not in PAR:
            raise Exception

        if 'DT' not in PAR:
            raise Exception

        if 'F0' not in PAR:
            raise Exception


    def check_solver_parameter_files(self):
        """ Checks solver parameters
        """
        nt = getpar('nt', cast=int)
        dt = getpar('deltat', cast=float)
        f0 = getpar('f0', file='DATA/SOURCE', cast=float)

        if nt != PAR.NT:
            if system.getnode() == 0: print "WARNING: nt != PAR.NT"
            setpar('nt', PAR.NT)

        if dt != PAR.DT:
            if system.getnode() == 0: print "WARNING: dt != PAR.DT"
            setpar('dt', PAR.DT)
        
        if f0 != PAR.F0:
            if system.getnode() == 0: print "WARNING: f0 != PAR.F0"
            setpar('f0', PAR.F0, file='DATA/SOURCE')

        if 'MULTIPLES' in PAR:
            if PAR.MULTIPLES:
                setpar('absorbtop', '.false.')
            else:
                setpar('absorbtop', '.true.')


    def generate_data(self, **model_kwargs):
        """ Generates data
        """
        self.generate_mesh(**model_kwargs)

        unix.cd(self.getpath)
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D')

        unix.mv(self.data_wildcard, 'traces/obs')
        self.export_traces(PATH.OUTPUT, 'traces/obs')


    def generate_mesh(self, model_path=None, model_name=None, model_type='gll'):
        """ Performs meshing and database generation
        """
        assert(model_name)
        assert(model_type)
        assert (exists(model_path))

        self.initialize_solver_directories()
        unix.cp(model_path, 'DATA/proc000000_rho_vp_vs.dat')
        self.export_model(PATH.OUTPUT +'/'+ model_name)


    def generate_precond(self, process_traces=None, model_path=None, model_name=None, model_type='gll'):
        assert(model_name)
        assert(model_type)
        assert (exists(model_path))

        self.initialize_solver_directories()
        unix.cp(model_path, 'DATA/proc000000_rho_vp_vs.dat')
        self.export_model(PATH.OUTPUT +'/'+ model_name)

        self.forward()
        unix.mv(self.data_wildcard, 'traces/syn')
        self.initialize_adjoint_traces('traces/syn')
        process_traces(self.getpath)

        self.adjoint()
        self.export_kernels(PATH.GLOBAL)


    ### low-level solver interface

    def forward(self):
        """ Calls SPECFEM2D forward solver
        """
        setpar('SIMULATION_TYPE', '1')
        setpar('SAVE_FORWARD', '.true.')
        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D')


    def adjoint(self):
        """ Calls SPECFEM2D adjoint solver
        """
        setpar('SIMULATION_TYPE', '3')
        setpar('SAVE_FORWARD', '.false.')
        unix.rm('SEM')
        unix.ln('traces/adj', 'SEM')

        self.mpirun('bin/xmeshfem2D')
        self.mpirun('bin/xspecfem2D')


    ### model input/output

    def load(self, filename, prefix='', suffix='', verbose=False):
        """Reads SPECFEM2D kernel or model

           Models and kernels are read from 5 or 6 column text files whose
           format is described in the SPECFEM2D user manual. Once read, a model
           or kernel is stored in a dictionary containing mesh coordinates and
           corresponding material parameter values.
        """
        # read text file
        M = np.loadtxt(filename)
        nrow = M.shape[0]
        ncol = M.shape[1]

        if ncol != 5:
            raise Exception('Bad SPECFEM2D model or kernel.')

        # fill in dictionary
        ioff = 0
        model = {}
        for key in ['x', 'z', 'rho', 'vp', 'vs']:
            model[key] = [M[:,ioff]]
            ioff += 1
        return model


    def save(self, filename, model, type='dummy', prefix='dummy', suffix='dummy'):
        """ writes SPECFEM2D kernel or model
        """
        nrow = len(model[model.keys().pop()][0])
        ncol = 5
        M = np.zeros((nrow, ncol))

        # fill in array
        for icol, key in enumerate(('x', 'z', 'rho', 'vp', 'vs')):
            if key in model.keys():
                M[:,icol] = model[key][0]
            else:
                M[:,icol] = loadbyproc(PATH.MODEL_INIT, key)

        # write array
        np.savetxt(filename, M, '%16.10e')



    ### postprocessing utilities

    def combine(self, path=''):
        """ Combines SPECFEM2D kernels
        """
        unix.cd(self.getpath)

        with open('kernel_paths', 'w') as f:
            f.writelines([join(path, dir)+'\n' for dir in unix.ls(path)])

        for name in self.parameters:
            self.mpirun(
                'bin/xcombine_rho_vp_vs '
                + 'kernel_paths' + ' '
                + path +'/'+ 'sum')


    def smooth(self, path='', span=0.):
        """ Smooths SPECFEM2D kernels by convolving them with a Gaussian
        """
        from seisflows.tools.array import meshsmooth, stack

        kernels = self.load(path, suffix='_kernel')
        if not span:
            return kernels

        x = kernels['x'][0]
        z = kernels['z'][0]
        mesh = stack(x, z)

        for key in self.parameters:
            kernels[key] = [meshsmooth(kernels[key][0], mesh, span)]

        unix.mv(path, path + '_nosmooth')
        self.save(path, kernels)


    def clip(self, path='', thresh=1.):
        """ Clips SPECFEM2D kernels
        """
        parts = self.load(path)
        if thresh >= 1.:
            return parts

        for key in self.parameters:
            # scale to [-1,1]
            minval = parts[key][0].min()
            maxval = parts[key][0].max()
            np.clip(parts[key][0], thresh*minval, thresh*maxval, out=parts[key][0])
        unix.mv(path, path + '_noclip')
        self.save(path, parts)


    ### file transfer utilities

    def import_model(self, path):
        src = join(path +'/'+ 'model')
        dst = join(self.getpath, 'DATA/proc000000_rho_vp_vs.dat')
        unix.cp(src, dst)

    def export_model(self, path):
        if system.getnode() == 0:
            src = join(self.getpath, 'DATA/proc000000_rho_vp_vs.dat')
            dst = path
            unix.cp(src, dst)

    def export_kernels(self, path):
        unix.mkdir_gpfs(join(path, 'kernels'))
        src = join(self.getpath, 'OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat')
        dst = join(path, 'kernels', '%06d' % system.getnode())
        unix.cp(src, dst)


    ### input file writers

    def write_parameters(self):
        unix.cd(self.getpath)
        solvertools.write_parameters(vars(PAR))

    def write_receivers(self):
        unix.cd(self.getpath)
        key = 'use_existing_STATIONS'
        val = '.true.'
        setpar(key, val)
        _, h = preprocess.load('traces/obs')
        solvertools.write_receivers(h.nr, h.rx, h.rz)

    def write_sources(self):
        unix.cd(self.getpath)
        _, h = preprocess.load(dir='traces/obs')
        solvertools.write_sources(vars(PAR), h)


    ### utility functions

    def mpirun(self, script, output='/dev/null'):
        """ Wrapper for mpirun
        """
        with open(output,'w') as f:
            subprocess.call(
                script,
                shell=True,
                stdout=f)

    ### miscellaneous

    @property
    def data_wildcard(self):
        return glob('OUTPUT_FILES/U?_file_single.su')

    @property
    def model_databases(self):
        return join(self.getpath, 'OUTPUT_FILES/DATABASES_MPI')

    @property
    def source_prefix(self):
        return 'SOURCE'


def loadbyproc(filename, key, nproc=None):
    # read text file
    M = np.loadtxt(filename)
    nrow = M.shape[0]
    ncol = M.shape[1]

    if ncol != 5:
        raise Exception('Bad SPECFEM2D model or kernel.')

    if key == 'x':
        return M[:, 0]
    elif key == 'z':
        return M[:, 1]
    elif key == 'rho':
        return M[:, 2]
    elif key == 'vp':
        return M[:, 3]
    elif key == 'vs':
        return M[:, 4]

