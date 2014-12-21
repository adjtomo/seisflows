import numpy as np

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.code import exists, glob, join
from seisflows.tools.config import loadclass, findpath, ConfigObj, ParameterObj
from seisflows.seistools.core import SeisStruct

OBJ = ConfigObj('SeisflowsObjects')
PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system
import preprocess


class specfem2d_SourceEncoding(loadclass('solver', 'specfem2d')):
    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(self.__class__, self).check()

        if 'NT_PADDED' not in PAR:
            raise Exception

    def prepare_dirs(self):
        """ Sets up directory in which to run solver
        """
        super(specfem2d_SourceEncoding, self).prepare_dirs()
        seistools.specfem2d.setpar('NSOURCES', PAR.NSRC)
        seistools.specfem2d.setpar('nt', PAR.NT_PADDED)

    def write_receivers(self):
        """ Writes receivers file
        """
        unix.cd(self.path)

        _, hdr = preprocess.load('traces/obs')
        seistools.specfem2d.write_receivers(hdr)

    def write_sources(self, sinfo=[], mapping=lambda i: [i]):
        """ Writes sources file
        """
        unix.cd(self.path)

        nodes = mapping(system.getnode())
        lines = []
        for i in nodes:
            seistools.specfem2d.write_sources(vars(PAR), sinfo[i])
            with open('DATA/SOURCE', 'r') as f:
                lines.extend(f.readlines())

        with open('DATA/SOURCE', 'w') as f:
            f.writelines(lines)

    def initialize_adjoint(self):
        zeros = np.zeros((PAR.NT_PADDED, PAR.NREC))
        h = SeisStruct(PAR.NREC, PAR.NT_PADDED, PAR.DT)
        for channel in ['x', 'y', 'z']:
            self.writer(zeros, h, channel=channel, prefix='traces/adj')
