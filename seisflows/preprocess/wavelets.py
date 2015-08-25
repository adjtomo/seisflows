
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class wavelets(loadclass('preprocess', 'base')):
    """ Data preprocessing class
    """
    raise NotImplementedError

