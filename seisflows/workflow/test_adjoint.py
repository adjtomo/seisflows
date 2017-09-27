
import sys
import numpy as np

from glob import glob
from os.path import basename, join
from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import ParameterError
from seisflows.workflow.base import base


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']
preprocess = sys.modules['seisflows_preprocess']


def DotProductLHS(keys, x, y):
    val = 0
    for key in keys:
        a = x[key].flatten()
        b = y[key].flatten()
        val += np.dot(a,b)
    val *= PAR.DT**2
    return val


def DotProductRHS(keys, x, y):
    val = 0
    for key in keys:
        a = np.array([])
        b = np.array([])
        for iproc in range(PAR.NPROC):
            a = np.append(a, x[key][iproc])
            b = np.append(b, y[key][iproc])
        val += np.dot(a,b)
    return val



class test_adjoint(base):

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

        if 'SOLVER' not in PATH:
            raise ParameterError(PATH, 'SOLVER')

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
        unix.rm(PATH.SCRATCH)
        unix.mkdir(PATH.SCRATCH)
        preprocess.setup()


        print 'SIMULATION 1 OF 3'
        system.run('solver', 'setup')

        print 'SIMULATION 2 OF 3'
        self.prepare_model()
        system.run('solver', 'eval_func',
                   path=PATH.SCRATCH)

        print 'SIMULATION 3 OF 3'
        system.run('solver', 'eval_grad',
                   path=PATH.SCRATCH)

        # collect traces
        obs = join(PATH.SOLVER, self.event, 'traces/obs')
        syn = join(PATH.SOLVER, self.event, 'traces/syn')
        adj = join(PATH.SOLVER, self.event, 'traces/adj')

        obs,_ = preprocess.load(obs)
        syn,_ = preprocess.load(syn)
        adj,_ = preprocess.load(adj, suffix='.su.adj')

        # collect model and kernels
        model = solver.load(PATH.MODEL_INIT)
        kernels = solver.load(PATH.SCRATCH+'/'+'kernels'+'/'+self.event, suffix='_kernel')

        # dot prodcut in data space
        keys = obs.keys()
        LHS = DotProductLHS(keys, syn, adj)

        # dot product in model space
        keys = ['rho', 'vp', 'vs'] # model.keys()
        RHS = DotProductRHS(keys, model, kernels)

        print 
        print 'LHS:', LHS
        print 'RHS:', RHS
        print 'RELATIVE DIFFERENCE:', (LHS-RHS)/RHS
        print


    ### utility functions

    def prepare_model(self):
        model = PATH.OUTPUT +'/'+ 'model_init'
        assert exists(model)
        unix.ln(model, PATH.SCRATCH +'/'+ 'model')

    @property
    def event(self):
        if not hasattr(self, '_event'):
            self._event = basename(glob(PATH.OUTPUT+'/'+'traces/*')[0])
        return self._event

