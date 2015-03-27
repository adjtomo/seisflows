
from glob import glob

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import solver
import preprocess
import postprocess


class migration(object):
    """ Migration base class.

      In the terminology of seismic exploration, implements a
      'reverse time migration'.
    """

    def check(self):
        """ Checks parameters and paths
        """
        # check paths
        if 'GLOBAL' not in PATH:
            raise ParameterError(PATH, 'GLOBAL')

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'OUTPUT' not in PATH:
            raise ParameterError(PATH, 'OUTPUT')

        # check input
        if 'DATA' not in PATH:
            setattr(PATH, 'DATA', None)

        if not exists(PATH.DATA):
            assert 'MODEL_TRUE' in PATH

        if 'MODEL_INIT' not in PATH:
            raise ParameterError(PATH, 'MODEL_INIT')

        # check output
        if 'SAVEGRADIENT' not in PAR:
            setattr(PAR, 'SAVEGRADIENT', 1)

        if 'SAVEKERNELS' not in PAR:
            setattr(PAR, 'SAVEKERNELS', 0)

        if 'SAVETRACES' not in PAR:
            setattr(PAR, 'SAVETRACES', 0)

        # assertions
        if 'OPTIMIZE' in PAR:
            assert not PAR.OPTIMIZE, "To run a migration, set PAR.OPTIMIZE = None"


    def main(self):
        """ Migrates seismic data
        """
        # prepare directory structure
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)

        # set up workflow machinery
        preprocess.setup()
        postprocess.setup()

        # set up solver machinery
        print 'Preparing solver...'
        system.run('solver', 'setup',
                   hosts='all')

        self.prepare_model()

        print 'Generating synthetics...'
        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.GLOBAL)

        # backproject data
        print 'Backprojecting data...'
        system.run('solver', 'eval_grad',
                   hosts='all',
                   path=PATH.GLOBAL,
                   export_traces=PAR.SAVETRACES)

        # process gradient
        postprocess.combine_kernels(
            path=PATH.GLOBAL)

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVEKERNELS:
            self.save_kernels()
        else:
            self.save_kernels_sum()

        print 'Finished\n'


    ### utility functions

    def prepare_model(self):
        model = PATH.OUTPUT +'/'+ 'model_init'
        assert exists(model)
        unix.ln(model, PATH.GLOBAL +'/'+ 'model')

    def save_kernels_sum(self):
        src = PATH.GLOBAL +'/'+ 'kernels/sum'
        dst = PATH.OUTPUT +'/'+ 'kernels/sum'
        unix.mkdir(dst)
        unix.mv(src, dst)

    def save_kernels(self):
        src = PATH.GLOBAL +'/'+ 'kernels'
        dst = PATH.OUTPUT
        unix.mkdir(dst)
        unix.mv(src, dst)

    def save_traces(self):
        src = PATH.GLOBAL +'/'+ 'traces'
        dst = PATH.OUTPUT
        unix.mv(src, dst)

