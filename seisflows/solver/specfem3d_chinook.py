#!/usr/bin/env python3
"""
This class provides utilities for the Seisflows solver interactions with
Specfem3D Cartesian, specifically to run on Chinook.
"""
import sys
import subprocess
from seisflows import logger
from seisflows.tools import unix, msg
from seisflows.solver.specfem3d import Specfem3D


class Specfem3D_Chinook(Specfem3D):
    """
    Solver SPECFEM3D
    ----------------
    SPECFEM3D-specific alterations to the base SPECFEM module

    Parameters
    ----------

    Paths
    -----
    ***
    """
    __doc__ = Specfem3D.__doc__ + __doc__

    def _run_binary(self, executable, stdout="solver.log", with_mpi=True):
        """
        Calls MPI solver executable to run solver binaries, used by individual
        processes to run the solver on system. If the external solver returns a
        non-zero exit code (failure), this function will return a negative
        boolean.

        .. note::
            This function ASSUMES it is being run from a SPECFEM working
            directory, i.e., that the executables are located in ./bin/

        .. note::
            This is essentially an error-catching wrapper of subprocess.run()

        :type executable: str
        :param executable: executable function to call. May or may not start
            E.g., acceptable calls for the solver would './bin/xspecfem2D'.
            Also accepts additional command line arguments such as:
            'xcombine_sem alpha_kernel kernel_paths...'
        :type stdout: str
        :param stdout: where to redirect stdout
        :type with_mpi: bool
        :param with_mpi: If `mpiexec` is given, use MPI to run the executable.
            Some executables (e.g., combine_vol_data_vtk) must be run in
            serial so this flag allows them to turn off MPI running.
        :raises SystemExit: If external numerical solver return any failure
            code while running
        """
        # Executable may come with additional sub arguments, we only need to
        # check that the actually executable exists
        if not unix.which(executable.split(" ")[0]):
            logger.critical(msg.cli(f"executable '{executable}' does not exist",
                            header="external solver error", border="="))
            sys.exit(-1)

        # Prepend with `mpiexec` if we are running with MPI
        # looks something like: `mpirun -n 4 ./bin/xspecfem2d`
        if self._mpiexec and with_mpi:
            executable = f"{self._mpiexec} -n {self.nproc} {executable}"
        logger.debug(f"running executable with cmd: '{executable}'")

        try:
            with open(stdout, "w") as f:
                subprocess.run(executable, shell=True, check=True, stdout=f,
                               stderr=f)
        except (subprocess.CalledProcessError, OSError) as e:
            logger.critical(
                msg.cli("The external numerical solver has returned a "
                        "nonzero exit code (failure). Consider stopping any "
                        "currently running jobs to avoid wasted "
                        "computational resources. Check 'scratch/solver/"
                        f"mainsolver/{stdout}' for the solvers stdout log "
                        "message. The failing command and error message are:",
                        items=[f"exc: {executable}", f"err: {e}"],
                        header="external solver error",
                        border="=")
            )
            sys.exit(-1)