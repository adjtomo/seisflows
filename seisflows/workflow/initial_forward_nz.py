import os
import sys
import time
import glob
import subprocess
import numpy as np

from seisflows.tools import unix
from seisflows.tools.tools import divides, exists
from seisflows.config import ParameterError, save
from seisflows.workflow.base import base

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']
optimize = sys.modules['seisflows_optimize']
preprocess = sys.modules['seisflows_preprocess']
postprocess = sys.modules['seisflows_postprocess']


class initial_forward_nz(base):
    """ 
    Forward run class

    A pared down version of inversion_nz

    Performs a suite of forward runs with misfit quantification.
    This is primarily useful for the very beginning of an inversion where
    manual review of misfit windows, data quality and event coverage is
    necessary. 
    """

    def check(self):
        """ Checks parameters and paths
        """

        # scratch paths
        if 'SCRATCH' not in PATH:
            raise ParameterError(PATH, 'SCRATCH')

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'FUNC' not in PATH:
            setattr(PATH, 'FUNC', os.path.join(PATH.SCRATCH, 'evalfunc'))

        if 'GRAD' not in PATH:
            setattr(PATH, 'GRAD', os.path.join(PATH.SCRATCH, 'evalgrad'))

        if 'OPTIMIZE' not in PATH:
            setattr(PATH, 'OPTIMIZE', os.path.join(PATH.SCRATCH, 'optimize'))

        # input paths
        if 'MODEL_INIT' not in PATH:
            raise ParameterError(PATH, 'MODEL_INIT')

        # output paths
        if 'OUTPUT' not in PATH:
            raise ParameterError(PATH, 'OUTPUT')

        if 'SAVEMODEL' not in PAR:
            setattr(PAR, 'SAVEMODEL', 1)

        if 'SAVETRACES' not in PAR:
            setattr(PAR, 'SAVETRACES', 0)

        # pyatoa specific paths
        # Config file should be present here.
        if 'PYATOA_IO' not in PATH:
            raise ParameterError(PATH, 'PYATOA_IO')

        # make sure the Pyatoa plugin run script is present
        if 'PYATOA_RUN' not in PATH:
            raise ParameterError(PATH, 'PYATOA_RUN')

        # make sure a Python3 binary is avilalable
        if 'PYTHON3' not in PATH:
            raise ParameterError(PATH, 'PYTHON3')

        if not exists(PATH.MODEL_INIT):
            raise Exception()

    def main(self):
        """ Carries out seismic inversion
        """
        self.setup()
        print ''
        
        self.initialize()
        self.finalize()
        self.clean()

    def setup(self):
        """ Lays groundwork for inversion
        """
        print 'SETUP'
        optimize.setup()

        print '\tInitializing Pyatoa'
        pyatoa_init = " ".join([
            PATH.PYTHON3,
            PATH.PYATOA_RUN,
            "--mode initialize",
            "--working_dir {}".format(PATH.WORKDIR)
        ])
        try:
            stdout = subprocess.check_output(pyatoa_init, shell=True)
        except subprocess.CalledProcessError as e:
            print("Pyatoa failed with {}".format(e))
            sys.exit(-1)

        print '\tPreparing initial model'
        system.run('solver', 'setup')

    def initialize(self):
        """ Prepares for next model update iteration
        """
        print 'INITIALIZE'
        self.write_model(path=PATH.GRAD, suffix='new')

        print '\tRunning forward simulation'
        print '\t\tstarting at', time.asctime(), '...'
        system.run('solver', 'eval_fwd', path=PATH.GRAD)
        print '\t\tfinished at', time.asctime()

        print '\tQuantifying misfit'
        print '\t\tstarting at', time.asctime(), '...'
        system.run_ancil('solver', 'eval_func',
                         iter=optimize.iter, suffix='new')
        print '\t\tfinished at', time.asctime()

        print '\tWriting misfit'
        self.write_misfit(suffix='new')

    def finalize(self):
        """ Saves results from current model update iteration
        """
        print 'FINALIZE'
        self.checkpoint()

        # Finalize Pyatoa run
        print '\tFinalizing Pyatoa'
        finalize_pyatoa = " ".join([
            PATH.PYTHON3,
            PATH.PYATOA_RUN,
            "--mode finalize",
            "--working_dir {}".format(PATH.WORKDIR),
            "--model_number {}".format("m{:0>2}".format(int(optimize.iter)-1)),
        ])
        try:
            stdout = subprocess.check_output(finalize_pyatoa, shell=True)
        except subprocess.CalledProcessError as e:
            print("Pyatoa failed with {}".format(e))
            sys.exit(-1)
            
        if divides(optimize.iter, PAR.SAVEMODEL):
            self.save_model()

        if divides(optimize.iter, PAR.SAVEGRADIENT):
            self.save_gradient()

        if divides(optimize.iter, PAR.SAVEKERNELS):
            self.save_kernels()

        if divides(optimize.iter, PAR.SAVETRACES):
            self.save_traces()

        if divides(optimize.iter, PAR.SAVERESIDUALS):
            self.save_residuals()

    def checkpoint(self):
        """ Writes information to disk so workflow can be resumed following a
          break
        """
        save()

    def write_model(self, path='', suffix=''):
        """ Writes model in format expected by solver
        Only write vp and vs, do not update density (rho) as it is not sensitive
        to the measurements we are making
        """
        src = 'm_'+suffix
        dst = path +'/'+ 'model'
        solver.save(dict=solver.split(optimize.load(src)), path=dst)

    def write_misfit(self, suffix=''):
        """ Writes misfit in format expected by nonlinear optimization library
            Overloads old write_misfit function
            Waits for all instances of Pyatoa to finish writing their misfit
        """
        src = os.path.join(PATH.PYATOA_IO, 'data', 'misfits', "*")
        dst = 'f_'+suffix
        
        misfit = 0
        while True:
            misfits = glob.glob(src)
            if len(misfits) == PAR.NSRC:
                for _misfit in misfits:
                    misfit += np.loadtxt(_misfit)
                    os.remove(_misfit)
                # Following Tape et al 2007, total misfit=misfit/number_srcs
                optimize.savetxt(dst, misfit/PAR.NSRC)
                return
            else:
                time.sleep(5)   

    def save_model(self):
        src = 'm_new'
        dst = os.path.join(PATH.OUTPUT, 'model_%04d' % optimize.iter)
        solver.save(solver.split(optimize.load(src)), dst)

    def save_traces(self):
        src = os.path.join(PATH.GRAD, 'traces')
        dst = os.path.join(PATH.OUTPUT, 'traces_%04d' % optimize.iter)
        unix.mv(src, dst)

