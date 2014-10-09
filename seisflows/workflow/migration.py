
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configtools import getclass, ParameterObject

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()
solver = getclass('solver',PAR.SOLVER)()



class migration(object):
    """ Migration base class.

      In the terminology of seismic exploration, implements a
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

        if 'BEGIN' not in PAR:
            setattr(PAR,'BEGIN',1)

        if 'END' not in PAR:
            setattr(PAR,'END',1)


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


    def main(self):
        """ Migrates seismic data
        """
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)


        # prepare solver
        print 'Preparing solver...'
        system.run( solver.prepare_solver,
            hosts='all' )

        self.prepare_model()

        system.run( solver.evaluate_func,
            hosts='all',
            path=PATH.SOLVER )


        # backproject data
        print 'Backprojecting data...'
        system.run( solver.evaluate_grad,
              hosts='all',
              path=PATH.SOLVER,
              export_traces=PAR.SAVETRACES )


        # process image
        self.postprocess = getclass('postprocess',PAR.POSTPROCESS)()
        self.postprocess.process_kernels(input=PATH.SOLVER)


        # save results
        if PAR.SAVEGRADIENT:
            self.save_gradient()

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVEKERNELS:
            self.save_kernels()

        print 'Finished\n'


    ### utility functions

    def prepare_model(self):
        model = PATH.OUTPUT+'/'+'model_init'
        assert exists(model)
        unix.cp(model,PATH.SOLVER+'/'+'model')

    def save_gradient(self):
        src = glob(join(PATH.POSTPROCESS,'gradient*'))
        dst = join(PATH.OUTPUT)
        unix.mv(src,dst)

    def save_kernels(self):
        src = join(PATH.SOLVER,'kernels')
        dst = join(PATH.OUTPUT)
        unix.mkdir(dst)
        unix.mv(src,dst)

    def save_traces(self):
        src = join(PATH.SOLVER,'traces')
        dst = join(PATH.OUTPUT)
        unix.mv(src,dst)

