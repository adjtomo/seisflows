
from glob import glob
from os.path import join

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import solver
import preprocess


class test_adjoint(object):

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


        # assertions
        if PAR.NSRC != 1:
            raise ParameterError(PAR, 'NSRC')


    def main(self):
        # prepare directory structure
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)

        # set up workflow machinery
        preprocess.setup()

        # set up solver machinery
        print 'Preparing solver...'
        system.run('solver', 'setup',
                   hosts='all')

        self.prepare_model()

        print 'Generating synthetics...'
        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.GLOBAL,
                   export_traces=True)

        print 'Backprojecting data...'
        system.run('solver', 'eval_grad',
                   hosts='all',
                   path=PATH.GLOBAL)

        # collect observations and synthetics
        src = join(PATH.OUTPUT, 'traces', '*')
        dst = join(PATH.OUTPUT, 'traces_obs')
        unix.mv(glob(src), dst)

        src = join(PATH.GLOBAL, 'traces', '*')
        dst = join(PATH.OUTPUT, 'traces_syn')
        unix.mv(glob(src), dst)

        obs = join(PATH.OUTPUT, 'traces_obs')
        syn = join(PATH.OUTPUT, 'traces_syn')

        d = preprocess.load(obs)
        s = preprocess.load(syn)

        model = solver.load(PATH.MODEL_INIT)
        kernels = solver.load(PATH.GLOBAL+'/'+'kernels/sum', suffix='_kernel')

        print 'd', d.keys()
        print 'd', d.keys()
        print 'model', model.keys()
        print 'kernels', kernels.keys()

        # LHS terms
        Ax = s
        y = d - s

        # RHS terms
        x = model
        Ay = kernels


    ### utility functions

    def prepare_model(self):
        model = PATH.OUTPUT +'/'+ 'model_init'
        assert exists(model)
        unix.ln(model, PATH.GLOBAL +'/'+ 'model')


