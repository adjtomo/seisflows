
from glob import glob

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

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
        if 'SCRATCH' not in PATH:
            raise ParameterError(PATH, 'SCRATCH')

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


    def main(self):
        """ Migrates seismic data
        """
        # prepare directory structure
        unix.rm(PATH.SCRATCH)
        unix.mkdir(PATH.SCRATCH)

        # set up workflow machinery
        preprocess.setup()
        postprocess.setup()

        # set up solver machinery
        print 'Preparing solver...'
        system.run('solver', 'setup',
                   hosts='all')

        self.prepare_model()

        # perform migration
        print 'Generating synthetics...'
        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.SCRATCH)

        print 'Backprojecting data...'
        system.run('solver', 'eval_grad',
                   hosts='all',
                   path=PATH.SCRATCH,
                   export_traces=PAR.SAVETRACES)

        postprocess.combine_kernels(
            path=PATH.SCRATCH,
            parameters=solver.parameters)

        try:
            postprocess.combine_kernels(
                path=PATH.SCRATCH,
                parameters=['rhop'])
        except:
            pass

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
        unix.cp(model, PATH.SCRATCH +'/'+ 'model')

    def save_kernels_sum(self):
        src = PATH.SCRATCH +'/'+ 'kernels/sum'
        dst = PATH.OUTPUT +'/'+ 'kernels'
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_kernels(self):
        src = PATH.SCRATCH +'/'+ 'kernels'
        dst = PATH.OUTPUT
        unix.mkdir(dst)
        unix.cp(src, dst)

    def save_traces(self):
        src = PATH.SCRATCH +'/'+ 'traces'
        dst = PATH.OUTPUT
        unix.cp(src, dst)

