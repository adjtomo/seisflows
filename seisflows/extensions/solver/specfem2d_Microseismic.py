import subprocess

import numpy as np
from scipy.interpolate import griddata

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join, setdiff, Struct
from seisflows.tools.configtools import getclass, getpath, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class specfem2d_Microseismic(getclass('solver','specfem2d')):

    # model parameters
    model_parameters = []
    model_parameters += ['rho']
    model_parameters += ['vp']
    model_parameters += ['vs']

    # inversion parameters
    inversion_parameters = []
    inversion_parameters += ['vs']

    # data channels
    channels = []
    channels += ['y']

    glob = lambda _ : glob('OUTPUT_FILES/*.semd')


    def __init__(self):
        super(specfem2d_Interferometry,self).__init__()
        self.preprocess.reader = seistools.specfem2d.read
        self.preprocess.writer = seistools.specfem2d.write


    def prepare_data(self,**kwargs):
        """ Prepares data for inversion or migration
        """
        #PAR.NT *= 2
        super(specfem2d_Interferometry,self).prepare_data(**kwargs)
        #PAR.NT /= 2


    def generate_data(self,**kwargs):
        """ Generates preprocess using SPECFEM2D forward solver
        """
        unix.cd(self.path)

        unix.mkdir('DATA/NOISE_TOMOGRAPHY')
        unix.mkdir('OUTPUT_FILES/NOISE_TOMOGRAPHY')

        # prepare model
        self.prepare_model(**kwargs)

        # copy input files
        src = glob(PATH.SOLVER_FILES+'/'+'*')
        dst = 'DATA/'
        unix.cp(src,dst)

        with open('DATA/NOISE_TOMOGRAPHY/irec_master','w') as f:
            f.write(str(system.getnode()+1)+'\n')

        # generating wavefield
        seistools.specfem2d.setpar('NOISE_TOMOGRAPHY', '1')
        seistools.specfem2d.setpar('SIMULATION_TYPE', '1')
        seistools.specfem2d.setpar('SAVE_FORWARD', '.false.')

        with open('log.fwd','w') as f:
            subprocess.call(self.mesher_binary,stdout=f)
            subprocess.call(self.solver_binary,stdout=f)

        # ensemble forward wavefield
        seistools.specfem2d.setpar('NOISE_TOMOGRAPHY', '2')
        seistools.specfem2d.setpar('SIMULATION_TYPE',  '1')
        seistools.specfem2d.setpar('SAVE_FORWARD', '.true.')

        with open('log.fwd','w') as f:
            subprocess.call(self.mesher_binary,stdout=f)
            subprocess.call(self.solver_binary,stdout=f)

        # export traces
        unix.mv(self.glob(),'traces/obs')


    ### high-level interface

    def evaluate_func(self,path='',export_traces=False):
        """ Evaluates misfit function by carrying out forward simulation and
            making measurements on observations and synthetics.
        """
        unix.cd(self.path)
        self.import_model(path)

        # forward simulation
        self.forward()
        unix.mv(self.glob(),'traces/syn')
        self.preprocess.prepare_adjoint(path=unix.pwd(),output_type=1)
        self.export_residuals(path)
        if export_traces:
            self.export_traces(path,'traces/syn')


    def evaluate_grad(self,path='',export_traces=False):
        """ Evaluates gradient by carrying out adjoint simulation. Adjoint traces
            must already be in place prior to calling this method or the adjoint
            simulation will fail.
        """
        unix.cd(self.path)

        # write residuals
        self.preprocess.prepare_adjoint(path=unix.pwd(),output_type=1)

        # negative branch
        self.prepare_adjoint(branch='neg')
        self.adjoint(tag='neg')

        # positive branch
        self.prepare_adjoint(branch='pos')
        self.adjoint(tag='pos')

        self.combine_pos_neg()
        self.export_kernels(path)



    def evaluate_hess(self,*args,**kwargs):
        raise NotImplementedError



    ### low-level interface

    def forward(self):
        """ Runs generating wavefield and ensembles forward wavefield
        """
        # adjust parameters
        seistools.specfem2d.setpar('NOISE_TOMOGRAPHY', '1')
        seistools.specfem2d.setpar('SIMULATION_TYPE',  '1')
        seistools.specfem2d.setpar('SAVE_FORWARD', '.false.')

        # generating wavefield
        with open('log.fwd','a') as f:
            subprocess.call(self.mesher_binary,stdout=f)
            subprocess.call(self.solver_binary,stdout=f)

        # adjust parameters
        seistools.specfem2d.setpar('NOISE_TOMOGRAPHY', '2')
        seistools.specfem2d.setpar('SIMULATION_TYPE',  '1')
        seistools.specfem2d.setpar('SAVE_FORWARD', ' .true.')

        # ensemble forward wavefield
        with open('log.fwd','w') as f:
            subprocess.call(self.mesher_binary,stdout=f)
            subprocess.call(self.solver_binary,stdout=f)


    def adjoint(self,tag='neg'):
        """ Calls adjoint solver
        """
        # adjust parameters
        seistools.specfem2d.setpar('NOISE_TOMOGRAPHY', '3')
        seistools.specfem2d.setpar('SIMULATION_TYPE',  '2')
        seistools.specfem2d.setpar('SAVE_FORWARD', '.false.')

        unix.rm('SEM')
        unix.ln('traces/adj','SEM')

        with open('log.fwd','w') as f:
            subprocess.call(self.mesher_binary,stdout=f)
            subprocess.call(self.solver_binary,stdout=f)

        src = 'OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat'
        dst = 'OUTPUT_FILES/kernel_'+tag
        unix.mv(src,dst)
        unix.rm(self.glob())



    def combine_pos_neg(self):
        """ Sums contributions from positive and negative branches
        """
        neg = self.load('OUTPUT_FILES/kernel_neg')
        pos = self.load('OUTPUT_FILES/kernel_pos')

        sum = {}
        sum['x'] = neg['x']
        sum['z'] = neg['z']

        for key in ['rho','vp','vs']:
            sum[key] = pos[key] + neg[key]

        file = 'OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat'
        self.save(file,sum,type='kernel')



    ### data processing utilities

    def prepare_adjoint(self,branch=''):
        d,h = self.preprocess.load(prefix='traces/obs')
        s,_ = self.preprocess.load(prefix='traces/syn')

        d = self.preprocess.apply(self.process_traces,[d],[h,branch])
        s = self.preprocess.apply(self.process_traces,[s],[h,branch])

        s = self.preprocess.apply(self.preprocess.compute_adjoint,[s,d],[h,2])
        self.preprocess.save(s,h)


    def process_traces(self,f,h,branch=''):
        """ Extracts positive or negative branch
        """
        i1 = 0
        i2 = int(np.ceil(h.nt/2.))

        i3 = int(np.floor(h.nt/2.))
        i4 = h.nt

        if branch == 'neg':
            f = seistools.swindow(f,h,i1,i2,units='samples')
            f[i1:i2,:] = np.flipud(f[i1:i2,:])
            f[i3:i4,:] = 0.
            f = self.preprocess.process_traces(f,h)
            f[i1:i2,:] = np.flipud(f[i1:i2,:])
            f[i3:i4,:] = 0.

        elif branch == 'pos':
            f = seistools.swindow(f,h,i3,i4,units='samples')
            f[i1:i2,:] = f[i3:i4,:]
            f[i3:i4,:] = 0.
            f = self.preprocess.process_traces(f,h)
            f[i3:i4,:] = f[i1:i2,:]
            f[i1:i2,:] = 0.

        return f



