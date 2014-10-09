
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import divides, exists, glob, irange, join
from seisflows.tools.configtools import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()



class GaussNewton(getclass('extensions.workflow','Newton')):
    """ Inversion with Gauss-Newton model updates
    """


    def __init__(self):
        """ Constructor
        """
        super(GaussNewton,self).__init__()


