
import sys
import numpy as np

from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.config import custom_import, ParameterError
from seisflows.workflow.base import base

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

solver = sys.modules['seisflows_solver']
postprocess = sys.modules['seisflows_postprocess']

migration = custom_import('workflow','migration')()


class test_postprocess(base):
    """ Postprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        migration.check()

        if 'INPUT' not in PATH:
            setattr(PATH, 'INPUT', None)


    def main(self):
        """ Writes gradient of objective function
        """
        if not PATH.INPUT:
            migration.main()

        postprocess.process_kernels()

