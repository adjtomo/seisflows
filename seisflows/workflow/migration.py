import numpy as np
import glob

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists, join
from seisflows.tools.config import loadclass, ParameterObj

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
        """ Checks parameters, paths, and dependencies
        """

        # check paths
        if 'GLOBAL' not in PATH:
            raise Exception

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'OUTPUT' not in PATH:
            raise Exception

        # check input
        if 'DATA' not in PATH:
            setattr(PATH, 'DATA', None)

        if not exists(PATH.DATA):
            assert 'MODEL_TRUE' in PATH

        if 'MODEL_INIT' not in PATH:
            raise Exception

        # check output
        if 'SAVEGRADIENT' not in PAR:
            setattr(PAR, 'SAVEGRADIENT', 1)

        if 'SAVEKERNELS' not in PAR:
            setattr(PAR, 'SAVEKERNELS', 0)

        if 'SAVETRACES' not in PAR:
            setattr(PAR, 'SAVETRACES', 0)


    def main(self):
        """ Migrates seismic data
        """
        # prepare directory structure
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)

        # set up pre- and post-processing
        preprocess.setup()
        postprocess.setup()

        # prepare solver
        print 'Preparing solver...'
        system.run('solver', 'setup',
                   hosts='all')

        self.prepare_model()

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
        postprocess.process_kernels(
            path=PATH.GLOBAL,
            tag='gradient')

        # save results
        if PAR.SAVEGRADIENT:
            self.save_gradient()

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVEKERNELS:
            self.save_kernels()

        print 'Finished\n'

    # -- utility functions

    def prepare_model(self):
        model = PATH.OUTPUT +'/'+ 'model_init'
        assert exists(model)
        unix.cp(model, PATH.GLOBAL +'/'+ 'model')

    def save_gradient(self):
        src = glob(PATH.GLOBAL +'/'+ 'gradient')
        dst = PATH.OUTPUT
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

