
import subprocess

from seisflows.tools import unix
from seisflows.tools.configtools import loadclass, ParameterObj

PAR = ParameterObj('parameters')
PATH = ParameterObj('paths')

system = loadclass('system',PAR.SYSTEM)()


class specfem3d_legacy(loadclass('solver','specfem3d')):

    def mpirun(self,script,outfile='/dev/null'):
        """ Wrapper for mpirun
        """
        unix.cd('bin')

        with open(outfile) as f:
            subprocess.call(
                  system.mpiargs()
                  + unix.basename(script),
                  shell=True,
                  stdout=f)
        unix.cd('..')
