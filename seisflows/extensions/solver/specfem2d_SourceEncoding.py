
import numpy as np

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configure import getclass, getpath, ParameterObject
from seisflows.seistools.core import SeisStruct

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class specfem2d_SourceEncoding(getclass('solver','specfem2d')):

    def write_parameters(self):
        """ Writes parameter file

            Calls base class method and then makes adjustments
        """
        super(specfem2d_SourceEncoding,self).write_parameters()

        seistools.specfem2d.setpar('NSOURCES',PAR.NSRC)
        seistools.specfem2d.setpar('nt',PAR.NT_PADDED)


    def write_sources(self,sinfo=[],mapping=lambda i:[i]):
        """ Writes sources file
        
            Writes sources file, taking into account that a single simulation 
            has to include multiple sources.
        """
        unix.cd(self.getpath())

        nodes = mapping(system.getnode())
        lines = []

        for i in nodes:
            seistools.specfem2d.write_sources(PAR.vars,sinfo[i])
            with open('DATA/SOURCE','r') as f:
                lines.extend(f.readlines())

        with open('DATA/SOURCE','w') as f:
            f.writelines(lines)


    def prepare_data(self,**kwargs):
        """ Prepares data for inversion or migration

            Overloads base method, removing write_sources and write_receivers
            since source and receiver factors are not yet available.
        """
        unix.cd(self.getpath())

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
            self.preprocess.writer(zeros,h,
                channel=channel,prefix='traces/adj')
