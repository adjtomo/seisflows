import os
import sys
import time
import glob
import subprocess
import numpy as np

from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import ParameterError, save, custom_import
from seisflows.workflow.base import base

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']
optimize = sys.modules['seisflows_optimize']
postprocess = sys.modules['seisflows_postprocess']


class until_adjoint(custom_import('workflow', 'inversion_nz')):
    """ 
    Run workflow for forward and adjoint simulations. Useful for creating
    waveform and kernel information, without invoking an entire inversion
    workflow. Supers 'inversion_nz'.
    """

    def check(self):
        """ Checks parameters and paths
        """
        # signifiy if data-synth. or synth.-synth. case
        if 'CASE' not in PAR:
            raise ParameterError(PAR, 'CASE')

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

        # check that their is a given starting model
        if not exists(PATH.MODEL_INIT):
            raise Exception()

        # synthetic-synthetic examples require a true model to create the 'data'
        if PAR.CASE == 'Synthetic' and not exists(PATH.MODEL_TRUE):
            raise Exception()

    def main(self):
        """ 
        Carries out forward and adjoint simulation
        """
        optimize.iter = 1
        main_start = time.time()
        print "BEGINNING WORKFLOW AT {}".format(time.asctime())
        self.setup()
        self.initialize()
        self.evaluate_gradient()

        print '{:.2f}m elapsed\n\n'.format((time.time() - main_start) / 60.)

    def setup(self):
        """ 
        Lays groundwork for inversion
        """
        postprocess.setup()
        optimize.setup()

        print 'SETUP'
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
        tstart = time.time()
        system.run('solver', 'setup')
        print '\t{:.2f}m elapsed'.format((time.time() - tstart) / 60.)


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
            "--step_count s00"
        ])
        try:
            stdout = subprocess.check_output(finalize_pyatoa, shell=True)
        except subprocess.CalledProcessError as e:
            print("Pyatoa failed with {}".format(e))
            sys.exit(-1)
