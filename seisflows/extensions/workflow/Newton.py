
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import divides, exists, glob, irange, join
from seisflows.tools.configtools import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()



class Newton(getclass('workflow','inversion')):
    """ Inversion with truncated Newton model updates
    """


    def __init__(self):
        """ Constructor
        """
        super(Newton,self).__init__()


    def compute_direction(self):
        """ Computes search direction
        """
        self.evaluate_gradient()

        self.optimize.initialize_newton()
        for self.ilcg in irange(1,PAR.LCGMAX):
            self.apply_hessian()
            isdone = self.optimize.update_newton()
            if isdone:
                break


    def apply_hessian(self):
        """ Computes the action of the Hessian on a given vector through
          appropriate solver call
        """
        if PAR.VERBOSE: print " LCG iteration", self.ilcg

        self.prepare_model(path=PATH.HESS, suffix='lcg')

        system.run( solver.apply_hess,
            hosts='all',
            path = PATH.HESS,
            hessian='Newton' )

        self.postprocess.process_kernels(
            input=PATH.HESS)

        src = PATH.POSTPROCESS+'/'+'gradient'
        g = solver.merge(solver.load(src))

        dst = PATH.OPTIMIZE+'/'+'g_lcg'
        savenpy(dst,g)

        unix.rm(PATH.HESS)
        unix.mkdir(PATH.HESS)

