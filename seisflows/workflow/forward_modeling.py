import numpy as np
import glob

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ConfigObj, ParameterObj

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class forward_modeling(object):
    """ Forward modeling base class
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """

        # check paths
        if 'GLOBAL' not in PATH:
            raise Exception

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        # check input settings
        if 'MODEL' not in PATH:
            raise Exception

        # check output settings
        if 'OUTPUT' not in PATH:
            raise Exception

        # check dependencies
        if 'solver' not in OBJ:
            raise Exception("Undefined Exception")

        if 'system' not in OBJ:
            raise Exception("Undefined Exception")

        global solver
        import solver

        global system
        import system

    def main(self):
        """ Generates seismic data
        """

        print 'Running solver...'

        system.run('solver', 'generate_data',
                   hosts='all',
                   model_path=PATH.MODEL,
                   model_type='gll',
                   model_name=None)

        print "Finished"
