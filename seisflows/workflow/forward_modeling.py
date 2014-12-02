
import numpy as np

from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configtools import loadclass, ConfigObj, ParameterObj

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
            setattr(PATH,'LOCAL',None)


        # check input settings
        if 'MODEL' not in PATH:
            raise Exception


        # check output settings
        if 'OUTPUT' not in PATH:
            raise Exception


        # check dependencies
        if 'solver' not in OBJ:
            raise Excpetion

        if 'system' not in OBJ:
            raise Excpetion

        global solver
        import solver

        global system
        import system



    def main(self):
        """ Generates seismic data
        """
       
        print 'Preparing solver directories...'

        system.run('solver','prepare_dirs',
          hosts = 'all')

        print 'Running solver...'

        system.run('solver','prepare_data',
          hosts = 'all',
          model_path = PATH.MODEL,
          model_type = 'gll',
          model_name = None)

        print "Finished"
