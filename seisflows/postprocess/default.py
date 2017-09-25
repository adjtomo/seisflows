
import sys

from seisflows.config import custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class default(custom_import('postprocess', 'base')):
    """ Default postprocesing option

      Provides default image processing and regularization functions for models
      or gradients
    """
    # currently identical to base class
    pass
