import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists, glob, join
from seisflows.tools.config import loadclass, ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


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

        if 'IMAGE' not in PATH:
            setattr(PATH, 'IMAGE', join(PATH.GLOBAL, 'image'))

        # check input paths
        if 'DATA' not in PATH:
            setattr(PATH, 'DATA', None)

        if not exists(PATH.DATA):
            assert 'MODEL_TRUE' in PATH

        if 'MODEL_INIT' not in PATH:
            raise Exception

        # check output paths
        if 'OUTPUT' not in PATH:
            raise Exception

        if 'SAVEIMAGE' not in PAR:
            setattr(PAR, 'SAVEIMAGE', 1)

        if 'SAVEKERNELS' not in PAR:
            setattr(PAR, 'SAVEKERNELS', 0)

        if 'SAVETRACES' not in PAR:
            setattr(PAR, 'SAVETRACES', 0)

        # check dependencies
        if 'postprocess' not in OBJ:
            raise Exception

        if 'solver' not in OBJ:
            raise Exception("Undefined Exception")

        if 'system' not in OBJ:
            raise Exception("Undefined Exception")

        global postprocess
        import postprocess

        global solver
        import solver

        global system
        import system

    def main(self):
        """ Migrates seismic data
        """
        # prepare directory structure
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)
        unix.mkdir(PATH.IMAGE)

        # prepare solver
        print 'Preparing solver...'
        system.run('solver', 'setup',
                   hosts='all')

        self.prepare_model()

        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.IMAGE)

        # backproject data
        print 'Backprojecting data...'
        system.run('solver', 'eval_grad',
                   hosts='all',
                   path=PATH.IMAGE,
                   export_traces=PAR.SAVETRACES)

        # process image
        postprocess.process_kernels(
            path=PATH.IMAGE,
            tag='image')

        # save results
        if PAR.SAVEIMAGE:
            self.save_image()

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVEKERNELS:
            self.save_kernels()

        print 'Finished\n'

    # -- utility functions

    def prepare_model(self):
        model = PATH.OUTPUT +'/'+ 'model_init'
        assert exists(model)
        unix.cp(model, PATH.IMAGE +'/'+ 'model')

    def save_image(self):
        src = glob(PATH.IMAGE +'/'+ 'image*')
        dst = PATH.OUTPUT
        unix.mv(src, dst)

    def save_kernels(self):
        src = PATH.IMAGE +'/'+ 'kernels'
        dst = PATH.OUTPUT
        unix.mkdir(dst)
        unix.mv(src, dst)

    def save_traces(self):
        src = PATH.IMAGE +'/'+ 'traces'
        dst = PATH.OUTPUT
        unix.mv(src, dst)

