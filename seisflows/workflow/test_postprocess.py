
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    loadclass, ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import solver
import postprocess

migration = loadclass('workflow','migration')()


class test_postprocess(object):
    """ Postprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        migration.check()

        if PATH.INPUT not in PAR:
            setattr(PATH, 'INPUT', None)


    def main(self):
        """ Writes gradient of objective function
        """
        if not PATH.INPUT:
            migration.main()

        postprocess.process_kernels()

