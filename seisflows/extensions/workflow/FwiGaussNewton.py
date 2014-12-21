import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import divides, exists, glob, irange, join
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('parameters')
PATH = ParameterObj('paths')


class FwiGaussNewton(loadclass('extensions.workflow', 'FwiNewton')):
    """ Inversion with Gauss-Newton model updates
    """

    raise NotImplementedError

