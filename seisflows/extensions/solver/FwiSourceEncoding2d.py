
import numpy as np

from seisflows import seistools
from seisflows.tools import unix
from seisflows.tools.codetools import exists, glob, join
from seisflows.tools.configure import getclass, getpath, ParameterObject
from seisflows.seistools.core import SeisStruct

PAR = ParameterObject('parameters')
PATH = ParameterObject('paths')

system = getclass('system',PAR.SYSTEM)()


class FwiSourceEncoding2d(getclass('solver','specfem2d')):

  def writerec(self):
      unix.cd(self.getpath())

      _,h = self.preprocess.load(prefix='traces/obs')
      seistools.specfem2d.writerec(h.nr,h.rx,h.rz)


  def writesrc(self,sinfo=[],mapping=lambda i:[i]):
      unix.cd(self.getpath())

      nodes = mapping(system.getnode())
      lines = []

      for i in nodes:
        seistools.specfem2d.writesrc(PAR.vars,sinfo[i])
        with open('DATA/SOURCE','r') as f:
          lines.extend(f.readlines())

      with open('DATA/SOURCE','w') as f:
        f.writelines(lines)


  def writepar(self):
      super(FwiSourceEncoding2d,self).writepar()

      seistools.specfem2d.setpar('NSOURCES',PAR.NSRC)
      seistools.specfem2d.setpar('nt',PAR.NT_PADDED)




  def prepare_data(self,**kwargs):
      """ Prepares data for inversion or migration

        Users implementing an inversion or migration can choose between 
        supplying data or supplying a target model from which data are 
        generated on the fly. If data are provided, SPECFEM2D input files will 
        be generated automatically from included header information. If data 
        are to be generated on the fly, both a target model and SPECFEM2D input
        files must be supplied.

        Adjoint traces are intialized by writing zeros for all components. This
        is necessary because, at the start of an adjoint simulation, SPECFEM2D
        expects all components to exist, even ones not actually in use for the
        inversion.
      """
      unix.cd(self.getpath())

      if exists(PATH.DATA):
          # generate SPECFEM2D input files
          self.writepar()

      else:
          # generate data
          self.generate_data(**kwargs)

      # initialize adjoint traces
      zeros = np.zeros((PAR.NT_PADDED,PAR.NREC))
      h = SeisStruct(PAR.NREC,PAR.NT_PADDED,PAR.DT)
      for channel in ['x','y','z']:
        self.preprocess.writer(zeros,h,
            channel=channel,prefix='traces/adj')


