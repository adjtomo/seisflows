
import numpy as np

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configtools import loadclass, findpath, ConfigObj, ParameterObj
from seisflows.seistools.core import SeisStruct

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class specfem2d_SourceEncoding(loadclass('solver','specfem2d')):

    def check(self):
        """ Checks objects and parameters
        """
        super(self.__class__,self).check()

        # check objects
        if 'system' not in OBJ:
            raise Excpetion

        global system
        import system

        # check parameters
        if 'NT_PADDED' not in PAR:
            raise Exception


    def write_parameters(self):
        """ Writes parameter file. Calls base method and makes adjustments
        """
        super(self.__class__,self).write_parameters()

        seistools.specfem2d.setpar('NSOURCES',PAR.NSRC)
        seistools.specfem2d.setpar('nt',PAR.NT_PADDED)


    def write_sources(self,sinfo=[],mapping=lambda i:[i]):
        """ Writes sources file
        """
        unix.cd(self.path)

        nodes = mapping(system.getnode())
        lines = []

        for i in nodes:
            seistools.specfem2d.write_sources(vars(PAR),sinfo[i])
            with open('DATA/SOURCE','r') as f:
                lines.extend(f.readlines())

        with open('DATA/SOURCE','w') as f:
            f.writelines(lines)


    def prepare_data(self,**kwargs):
        """ Postpones write_sources and write_receivers because required
            encoding factors are not yet available
        """
        unix.cd(self.path)

        if exists(PATH.DATA):
            # generate SPECFEM2D input files
            self.write_parameters()

        else:
            # generate data
            self.generate_data(**kwargs)

        # initialize adjoint traces
        zeros = np.zeros((PAR.NT_PADDED,PAR.NREC))
        h = SeisStruct(PAR.NREC,PAR.NT_PADDED,PAR.DT)
        for channel in ['x','y','z']:
            self.writer(zeros,h,
                channel=channel,prefix='traces/adj')
