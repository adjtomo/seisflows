import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import exists, glob, join
from seisflows.tools.config import loadclass, findpath, ParameterObj
from seisflows.seistools.shared import SeisStruct

import seisflows.seistools.specfem2d as solvertools

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

    def initialize_solver_directories(self):
        """ Sets up directory in which to run solver
        """
        super(specfem2d_SourceEncoding, self).initialize_solver_directories()
        solvertools.setpar('NSOURCES', PAR.NSRC)
        solvertools.setpar('nt', PAR.NT_PADDED)

    def write_receivers(self):
        """ Writes receivers file
        """
        unix.cd(self.spath)

        _, hdr = preprocess.load('traces/obs')
        solvertools.write_receivers(hdr)

    def write_sources(self, sinfo=[], mapping=lambda i: [i]):
        """ Writes sources file
        """
        unix.cd(self.spath)

        nodes = mapping(system.getnode())
        lines = []
        for i in nodes:
            solvertools.write_sources(vars(PAR), sinfo[i])
            with open('DATA/SOURCE', 'r') as f:
                lines.extend(f.readlines())

        with open('DATA/SOURCE', 'w') as f:
            f.writelines(lines)

    def initialize_adjoint_traces(self):
        zeros = np.zeros((PAR.NT_PADDED, PAR.NREC))
        h = SeisStruct(PAR.NREC, PAR.NT_PADDED, PAR.DT)
        for channel in ['x', 'y', 'z']:
            preprocess.writer(zeros, h, channel=channel, prefix='traces/adj')
