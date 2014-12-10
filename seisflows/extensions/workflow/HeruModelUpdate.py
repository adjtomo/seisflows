
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import divides, exists, glob, irange, join
from seisflows.tools.config import loadclass, ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class HeruModelUpdate(loadclass('workflow','inversion')):
    """ Given search direction, computes model update
    """

    def check(self):
        """ Checks parameters, paths and other prerequisites
        """

        # check dependencies
        if 'optimize' not in OBJ:
            raise Exception

        if 'postprocess' not in OBJ:
            raise Exception

        if 'solver' not in OBJ:
            raise Excpetion

        if 'system' not in OBJ:
            raise Excpetion

        global optimize
        import optimize

        global postprocess
        import postprocess

        global solver
        import solver

        global system
        import system


        # check user supplied parameters
        if PAR.BEGIN != 1:
            raise Exception

        if PAR.END != 1:
            raise Exception

        if PAR.SCHEME != 'GradientDescent':
            raise Exception

        if 'VERBOSE' not in PAR:
            setattr(PAR,'VERBOSE',1)


        # check scratch paths
        if 'GLOBAL' not in PATH:
            raise Exception

        if 'LOCAL' not in PATH:
            setattr(PATH,'LOCAL',None)

        if 'FUNC' not in PATH:
            setattr(PATH,'FUNC',join(PATH.GLOBAL,'func'))

        if 'GRAD' not in PATH:
            setattr(PATH,'GRAD',join(PATH.GLOBAL,'grad'))

        if 'HESS' not in PATH:
            setattr(PATH,'HESS',join(PATH.GLOBAL,'hess'))

        if 'OPTIMIZE' not in PATH:
            setattr(PATH,'OPTIMIZE',join(PATH.GLOBAL,'optimize'))


        # check input paths
        if 'DATA' not in PATH:
            setattr(PATH,'DATA',None)

        if not exists(PATH.DATA):
            assert 'MODEL_TRUE' in PATH

        if 'MODEL_INIT' not in PATH:
            raise Exception


        # check output paths
        if 'OUTPUT' not in PATH:
            raise Exception

        if 'SAVEMODEL' not in PAR:
            setattr(PAR,'SAVEMODEL',1)

        if 'SAVEGRADIENT' not in PAR:
            setattr(PAR,'SAVEGRADIENT',0)

        if 'SAVEKERNELS' not in PAR:
            setattr(PAR,'SAVEKERNELS',0)

        if 'SAVETRACES' not in PAR:
            setattr(PAR,'SAVETRACES',0)

        if 'SAVERESIDUALS' not in PAR:
            setattr(PAR,'SAVERESIDUALS',0)


        # check model update prerequisites
        assert not PATH.LOCAL
        assert exists(PATH.GRAD+'/'+'kernels/sum')
        assert exists(PATH.OPTIMIZE+'/'+'f_new')

        super(HeruModelUpdate,self).check()



    def main(self):
        self.setup()

        print "Computing search direction"
        self.compute_direction()

        print "Computing step length"
        self.line_search()

        self.finalize()
        print ''



    def setup(self):
        self.iter = 1
        optimize.iter = 1

        optimize.setup()



    def compute_direction(self):
        """ Given gradient, computes search direction
        """
        system.run( 'postprocess','process_kernels',
            hosts='head',
            path=PATH.GRAD,
            optim_path=PATH.OPTIMIZE+'/'+'g_new')

        optimize.compute_direction()



    def finalize(self):
        """ Saves results from most recent model update iteration
        """
        if divides(self.iter,PAR.SAVEMODEL):
            self.save_model()

        if divides(self.iter,PAR.SAVEKERNELS):
            self.save_kernels()

        if divides(self.iter,PAR.SAVETRACES):
            self.save_traces()
