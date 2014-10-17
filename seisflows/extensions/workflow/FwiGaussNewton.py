
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import divides, exists, glob, irange, join
from seisflows.tools.configtools import getclass, GlobalStruct

PAR = GlobalStruct('parameters')
PATH = GlobalStruct('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()



class FwiGaussNewton(getclass('extensions.workflow','FwiNewton')):
    """ Inversion with Gauss-Newton model updates
    """


    def __init__(self):
        """ Constructor
        """
        super(FwiGaussNewton,self).__init__()


