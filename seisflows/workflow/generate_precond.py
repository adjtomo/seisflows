
from glob import glob

import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import ParameterObj

from seisflows.seistools import adjoint

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import solver
import preprocess
import postprocess


class generate_precond(object):

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

        if 'MODEL_INIT' not in PATH:
            raise ParameterError(PATH, 'MODEL_INIT')


    def main(self):
        # prepare directory structure
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)

        # set up workflow machinery
        preprocess.setup()
        postprocess.setup()

        print 'Generating preconditioner...'
        system.run('solver', 'generate_precond',
                   hosts='all',
                   process_traces=process_traces,
                   model_path=PATH.MODEL_INIT,
                   model_name='model',
                   model_type='gll')

        postprocess.combine_kernels(
            path=PATH.GLOBAL)

        self.process_kernels(
            path=PATH.GLOBAL)

        # save preconditioner
        src = PATH.GLOBAL +'/'+ 'kernels/absval'
        dst = PATH.OUTPUT +'/'+ 'precond'
        unix.cp(src, dst)

        print 'Finished\n'


    def process_kernels(self, path):
        assert (exists(path))

        # take absolute value
        parts = solver.load(path +'/'+ 'kernels/sum')
        for key in solver.parameters:
            parts[key] = np.abs(parts[key])

        solver.save(path +'/'+ 'kernels/absval_noclip',
                    parts,
                    suffix='_kernel')

        # clip
        for key in solver.parameters:
            maxval = np.max(parts[key])
            if maxval > 0:
                parts[key] = np.clip(parts[key]/maxval, 0., 0.0333)

        solver.save(path +'/'+ 'kernels/absval',
                    parts,
                    suffix='_kernel')

        # smooth
        if PAR.SMOOTH > 0.:
            system.run('solver', 'smooth',
                       hosts='head',
                       path=path +'/'+ 'kernels/absval',
                       span=PAR.SMOOTH)


def process_traces(path):
    unix.cd(path)

    s, h = preprocess.load(prefix='traces/syn/')
    s = preprocess.apply(preprocess.process_traces, [s], [h])

    s = preprocess.apply(adjoint.prepare_precond, [s], [h])
    preprocess.save(s, h, prefix='traces/adj/')

