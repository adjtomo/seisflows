
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import loadclass, ParameterObj
from seisflows.optimize import lib

PAR = ParameterObj('parameters')
PATH = ParameterObj('paths')



class Newton(loadclass('optimize','default')):
    """ Implements truncated Newton algorithm
    """

    def __init__(cls):
        super(Newton,cls).__init__()

        if PAR.SCHEME != 'Newton':
            setattr(PAR,'SCHEME','Newton')

        if 'LCGMAX' not in PAR:
            setattr(PAR,'LCGMAX',2)

        if 'LCGTHRESH' not in PAR:
            setattr(PAR,'LCGTHRESH',np.inf)


    def setup(cls):
        super(Newton,cls).setup()
        cls.LCG = lib.LCG(cls.path,PAR.LCGTHRESH,PAR.LCGMAX)


    def initialize_newton(cls):
        """ description goes here
        """
        unix.cd(cls.path)
        cls.iter += 1

        m = loadnpy('m_new')
        p = cls.LCG.initialize()
        cls.delta = 1.e-3/max(abs(p))
        savenpy('m_lcg',m+p*cls.delta)


    def update_newton(cls):
        """ description goes here
        """
        unix.cd(cls.path)

        g0 = loadnpy('g_new')
        g = loadnpy('g_lcg')

        p,isdone = cls.LCG.update((g-g0)/cls.delta)

        if isdone:
            savenpy('p_new',p)
            savetxt('s_new',cls.delta*np.dot(g0,p))
        else:
            m = loadnpy('m_new')
            cls.delta = 1.e-3/max(abs(p))
            savenpy('m_lcg',m+p*cls.delta)

        return isdone

