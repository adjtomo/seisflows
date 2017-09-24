
import sys

from seisflows.config import custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class default(custom_import('preprocess', 'base')):
    """ Default preprocesing class

      Provides data processing functions for seismic traces, with options for
      data misfit, filtering, normalization and muting
    """
    # currently identical to base class
    pass

