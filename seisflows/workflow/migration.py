
import numpy as np

from seisflows.tools import unix
from seisflows.tools.arraytools import loadnpy, savenpy
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configtools import loadclass, ParameterObj

PAR = ParameterObj('parameters')
PATH = ParameterObj('paths')

system = loadclass('system',PAR.SYSTEM)()
solver = loadclass('solver',PAR.SOLVER)()



class migration(object):
    """ Migration base class.

      In the terminology of seismic exploration, implements a
      'reverse time migration'.
    """

    def __init__(self):

        # check scratch paths
        if 'GLOBAL' not in PATH:
            raise Exception

        if 'LOCAL' not in PATH:
            setattr(PATH,'LOCAL',None)

        if 'IMAGE' not in PATH:
            setattr(PATH,'IMAGE',join(PATH.GLOBAL,'image'))


        # check input paths
        if 'DATA' not in PATH:
            setattr(PATH,'DATA',None)

        if not exists(PATH.DATA):
            assert 'MODEL_TRUE' in PATH

        if 'MODEL_INIT' not in PATH:
            raise Exception


        # check output paths
        if 'OUTPUT' not in PATH:
            setattr(PATH,'OUTPUT',join(PATH.SUBMIT,'output'))

        if 'SAVEIMAGE' not in PAR:
            setattr(PAR,'SAVEIMAGE',1)

        if 'SAVEKERNELS' not in PAR:
            setattr(PAR,'SAVEKERNELS',0)

        if 'SAVETRACES' not in PAR:
            setattr(PAR,'SAVETRACES',0)


    def main(self):
        """ Migrates seismic data
        """
        # prepare directory structure
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)
        unix.mkdir(PATH.IMAGE)
        unix.mkdir(PATH.OUTPUT)

        # prepare solver
        print 'Preparing solver...'
        system.run( solver.prepare_solver,
            hosts='all' )

        self.prepare_model()

        system.run( solver.evaluate_func,
            hosts='all',
            path=PATH.IMAGE )

        # backproject data
        print 'Backprojecting data...'
        system.run( solver.evaluate_grad,
              hosts='all',
              path=PATH.IMAGE,
              export_traces=PAR.SAVETRACES )

        # process image
        self.postprocess = loadclass('postprocess',PAR.POSTPROCESS)()

        self.postprocess.process_kernels(
            path=PATH.IMAGE,
            tag='image')


        # save results
        if PAR.SAVEIMAGE:
            self.save_image()

        if PAR.SAVETRACES:
            self.save_traces()

        if PAR.SAVEKERNELS:
            self.save_kernels()

        print 'Finished\n'


    ### utility functions

    def prepare_model(self):
        model = PATH.OUTPUT+'/'+'model_init'
        assert exists(model)
        unix.cp(model,PATH.IMAGE+'/'+'model')

    def save_image(self):
        src = glob(join(PATH.IMAGE,'image*'))
        dst = join(PATH.OUTPUT)
        unix.mv(src,dst)

    def save_kernels(self):
        src = join(PATH.IMAGE,'kernels')
        dst = join(PATH.OUTPUT)
        unix.mkdir(dst)
        unix.mv(src,dst)

    def save_traces(self):
        src = join(PATH.IMAGE,'traces')
        dst = join(PATH.OUTPUT)
        unix.mv(src,dst)

