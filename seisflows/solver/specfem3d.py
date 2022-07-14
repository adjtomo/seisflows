#!/usr/bin/env python3
"""
This is the subclass seisflows.solver.Specfem3D
This class provides utilities for the Seisflows solver interactions with
Specfem3D Cartesian. It inherits all attributes from seisflows.solver.Base,
and overwrites these functions to provide specified interaction with Specfem3D
"""
import os

from seisflows.tools import unix
from seisflows.tools.specfem import setpar, getpar
from seisflows.solver.specfem import Specfem


class Specfem3D(Specfem):
    """
    SPECFEM3D_Cartesian specific parameters

    :type source_prefix: str
    :param source_prefix: Prefix of source files in path SPECFEM_DATA. Defaults
        to 'CMTSOLUTION'
    :type multiples: bool
    :param multiples: set an absorbing top-boundary condition
    """
    __doc__ = Specfem.__doc__ + __doc__

    def __init__(self, source_prefix="CMTSOLUTION", **kwargs):
        """Instantiate a Specfem3D_Cartesian solver interface"""

        super().__init__(source_prefix=source_prefix, **kwargs)

        # Define parameters based on material type
        if self.materials.upper() == "ACOUSTIC":
            self._parameters.append("vp")
        elif self.materials.upper() == "ELASTIC":
            self._parameters.append("vp")
            self._parameters.append("vs")

        # Overwriting the base class parameters
        self._acceptable_source_prefixes = ["CMTSOLUTION", "FORCESOLUTION"]
        self._required_binaries = ["xspecfem3D", "xmeshfem3D",
                                   "xgenerate_databases" "xcombine_sem",
                                   "xsmooth_sem"]

    def data_wildcard(self, comp="?"):
        """
        Returns a wildcard identifier for synthetic data

        TODO where does SU put its component?

        :rtype: str
        :return: wildcard identifier for channels
        """
        if self.data_format.upper() == "SU":
            return f"*_d?_SU"
        elif self.data_format.upper() == "ASCII":
            return f"*.?X{comp}.sem?"

    @property
    def model_databases(self):
        """
        The location of databases for model outputs, usually
        OUTPUT_FILES/DATABASES_MPI. Value is grabbed from the Par_file
        """
        local_path = getpar(key="LOCAL_PATH",
                            file=os.path.join(self.cwd, "DATA", "Par_file"))[1]
        return local_path

    @property
    def kernel_databases(self):
        """
        The location of databases for kernel outputs, usually the same as
        'model_databases'
        """
        return self.model_databases

    def forward_simulation(self, executables=None, save_traces=False,
                           export_traces=False):
        """
        Calls SPECFEM3D forward solver, exports solver outputs to traces dir

        :type executables: list or None
        :param executables: list of SPECFEM executables to run, in order, to
            complete a forward simulation. This can be left None in most cases,
            which will select default values based on the specific solver
            being called (2D/3D/3D_GLOBE). It is made an optional parameter
            to keep the function more general for inheritance purposes.
        :type save_traces: str
        :param save_traces: move files from their native SPECFEM output location
            to another directory. This is used to move output waveforms to
            'traces/obs' or 'traces/syn' so that SeisFlows knows where to look
            for them, and so that SPECFEM doesn't overwrite existing files
            during subsequent forward simulations
        :type export_traces: str
        :param export_traces: export traces from the scratch directory to a more
            permanent storage location. i.e., copy files from their original
            location
        """
        if executables is None:
            executables = ["bin/xgenerate_databases", "bin/xspecfem3D"]

        unix.cd(self.cwd)

        # SPECFEM3D has to deal with attenuation
        if self.attenuation:
            setpar(key="ATTENUATION", val=".true.", file="DATA/Par_file")
        else:
            setpar(key="ATTENUATION", val=".false`.", file="DATA/Par_file")

        super().forward_simulation(executables=executables,
                                   save_traces=save_traces,
                                   export_traces=export_traces
                                   )

    def adjoint_simulation(self, executables=None, save_kernels=False,
                           export_kernels=False):
        """
        Calls SPECFEM3D adjoint solver, creates the `SEM` folder with adjoint
        traces which is required by the adjoint solver

        :type executables: list or None
        :param executables: list of SPECFEM executables to run, in order, to
            complete an adjoint simulation. This can be left None in most cases,
            which will select default values based on the specific solver
            being called (2D/3D/3D_GLOBE). It is made an optional parameter
            to keep the function more general for inheritance purposes.
        :type save_kernels: str
        :param save_kernels: move the kernels from their native SPECFEM output
            location to another path. This is used to move kernels to another
            SeisFlows scratch directory so that they are discoverable by
            other modules. The typical location they are moved to is
            path_eval_grad
        :type export_kernels: str
        :param export_kernels: export/copy/save kernels from the scratch
            directory to a more permanent storage location. i.e., copy files
            from their original location. Note that kernel file sizes are LARGE,
            so exporting kernels can lead to massive storage requirements.
        """
        if executables is None:
            executables = ["bin/xspecfem3D"]

        # Make sure attenuation is OFF, if ON you'll get a floating point error
        unix.cd(self.cwd)
        setpar(key="ATTENUATION", val=".false.", file="DATA/Par_file")

        super().adjoint_simulation(executables=executables,
                                   save_kernels=save_kernels,
                                   export_kernels=export_kernels)
