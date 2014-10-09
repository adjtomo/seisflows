
import numpy as np

from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configtools import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()


class modeling(object):
    """ Forward modeling base class
    """

    def __init__(self):
        """ Constructor
        """
        # check user supplied paths
        if 'MODEL' not in PATH:
            setattr(PATH,'MODEL','')
            print 'Warning: PATH.MODEL not defined.'

        # configure parameters
        PAR.OPTIMIZE = None
        PAR.POSTPROCESS = None

        # configure paths
        PATH.OUTPUT = join(PATH.SUBMIT_DIR,'output')
        unix.mkdir(PATH.OUTPUT)

        PATH.SCRATCH = join(PATH.GLOBAL,'scratch')
        if PATH.LOCAL: PATH.SCRATCH = join(PATH.LOCAL,'scatch')


    def main(self):
        """ Generates seismic data
        """
        system.run( solver.prepare_solver,
          hosts = 'all',
          inversion = False)

        print "Generating data"
        system.run( solver.generate_data,
          hosts = 'all',
          model_path = PATH.MODEL,
          model_type = 'gll',
          model_name = None)

        print "Finished"
