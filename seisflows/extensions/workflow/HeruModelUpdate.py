
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import divides, exists, glob, irange, join
from seisflows.tools.configtools import getclass, GlobalStruct

PAR = GlobalStruct('parameters')
PATH = GlobalStruct('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()



class HeruModelUpdate(getclass('workflow','inversion')):
    """ Given search direction, computes model update
    """


    def __init__(self):
        """ Checks parameters, paths and other prerequisites
        """
        self.iter = 1

        # check user supplied parameters
        if PAR.BEGIN != 1:
            raise Exception

        if PAR.END != 1:
            raise Exception

        if PAR.SCHEME != 'GradientDescent':
            raise Exception

        if 'SAVEMODELS' not in PAR:
            setattr(PAR,'SAVEMODELS',1)

        if 'SAVEKERNELS' not in PAR:
            setattr(PAR,'SAVEKERNELS',0)

        if 'SAVETRACES' not in PAR:
            setattr(PAR,'SAVETRACES',0)


        # add some additional paths
        PATH.OUTPUT = join(PATH.SUBMIT_DIR,'output')
        unix.mkdir(PATH.OUTPUT)

        PATH.SCRATCH = join(PATH.GLOBAL,'scratch')
        if PATH.LOCAL: PATH.SCRATCH = join(PATH.LOCAL,'scatch')

        PATH.FUNC = join(PATH.SOLVER,'func')
        PATH.GRAD = join(PATH.SOLVER,'grad')
        PATH.HESS = join(PATH.SOLVER,'hess')


        # check model update prerequisites
        assert exists(PATH.GRAD+'/'+'kernels/sum')
        assert exists(PATH.OPTIMIZE+'/'+'f_new')



    def main(self):
        self.setup()

        print "Computing search direction"
        self.compute_direction()

        print "Computing step length"
        self.line_search()

        self.finalize()
        print ''



    def setup(self):
        super(HeruModelUpdate,self).setup()



    def compute_direction(self):
        """ Given gradient, computes search direction
        """
        self.postprocess.process_kernels()
        self.optimize.compute_direction()



    def finalize(self):
        """ Saves results from most recent model update iteration
        """
        if divides(self.iter,PAR.SAVEMODELS):
            self.save_model()

        if divides(self.iter,PAR.SAVEKERNELS):
            self.save_kernels()

        if divides(self.iter,PAR.SAVETRACES):
            self.save_traces()
