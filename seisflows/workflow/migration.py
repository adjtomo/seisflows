
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configure import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()



class migration(getclass('workflow','inversion')):
    """ Migration base class.

      Implements migration by subclassing inversion and overloading 'main'.

      In the terminology of seismic exploration, implements
      'reverse time migration'.
    """

    def __init__(self):
        """ Class constructor
        """
        self.iter = 0

        # check user supplied parameters
        if 'SAVEGRADIENT' not in PAR:
            setattr(PAR,'SAVEGRADIENT',1)

        if 'SAVETRACES' not in PAR:
            setattr(PAR,'SAVETRACES',0)

        if 'SAVEKERNELS' not in PAR:
            setattr(PAR,'SAVEKERNELS',0)

        # check user supplied paths
        if 'DATA' not in PATH:
            setattr(PATH,'DATA','')

        if not exists(PATH.DATA):
            assert 'MODEL_TRUE' in PATH

        if 'MODEL_INIT' not in PATH:
            raise Exception

        # add paths to global dictionary
        PATH.OUTPUT = join(PATH.SUBMIT_DIR,'output')
        unix.mkdir(PATH.OUTPUT)

        PATH.SCRATCH = join(PATH.GLOBAL,'scratch')
        if PATH.LOCAL: PATH.SCRATCH = join(PATH.LOCAL,'scatch')

        PATH.FUNC = join(PATH.SOLVER,'func')
        PATH.GRAD = join(PATH.SOLVER,'grad')
        PATH.HESS = join(PATH.SOLVER,'hess')


    def main(self):
        """ Migrates seismic data
        """
        self.setup()

        self.initialize()

        print "Computing search direction"
        self.compute_direction()

        self.finalize()
        print 'Finished\n'


    def finalize(self):
        """ Saves results from migration
        """
        if PAR.SAVEGRADIENT:
            self.save_gradient()

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVEKERNELS:
            self.save_kernels()


