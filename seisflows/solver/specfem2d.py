#!/usr/bin/env python3
"""
This is the subclass seisflows.solver.specfem2d

This class provides utilities for the Seisflows solver interactions with
Specfem2D. It inherits all attributes from seisflows.solver.Base,
"""
import os
from glob import glob

from seisflows.solver.specfem import Specfem
from seisflows.tools import unix
from seisflows.tools.specfem import getpar, setpar


class Specfem2D(Specfem):
    """
    Python interface to Specfem2D.
    """
    def __init__(self, source_prefix=None, multiples=False, **kwargs):
        """
        SPECFEM2D specific parameters

        :type source_prefix: str
        :param source_prefix: Prefix of SOURCE files in path SPECFEM_DATA.
        :type multiples: bool
        :param multiples: set an absorbing top-boundary condition
        """
        super().__init__(**kwargs)

        self.source_prefix = source_prefix or "SOURCE"
        self.multiples = multiples
        self.f0 = None

        # Define parameters based on material type
        if self.materials.upper() == "ACOUSTIC":
            self._parameters.append("vp")
        elif self.materials.upper() == "ELASTIC":
            self._parameters.append("vp")
            self._parameters.append("vs")

        self._acceptable_source_prefix = ["SOURCE", "FORCE", "FORCESOLUTION"]

    def setup(self):
        """
        Additional SPECFEM2D setup steps
        """
        super().setup()

        self.f0 = getpar(key="f0",
                         file=os.path.join(self.cwd, "DATA/SOURCE"))[1]

        if self.multiples:
            setpar(key="absorbtop", val=".false.", file="DATA/Par_file")
        else:
            setpar(key="absorbtop", val=".true.", file="DATA/Par_file")

    def smooth(self, input_path, output_path, parameters=None, span_h=0.,
               span_v=0., use_gpu=False, output="smooth.log"):
        """
        Specfem2D requires additional model parameters in directory to perform
        the xsmooth_sem task. This function will copy these files into the
        directory before performing the base smooth operations.

        Kwargs should match arguments of solver.base.smooth()

        .. note::
            This operation is usually run with run(single=True) so only one
            task will be performing these operations.

        :type input_path: str
        :param input_path: path to data
        :type output_path: str
        :param output_path: path to export the outputs of xcombine_sem
        :type parameters: list
        :param parameters: optional list of parameters,
            defaults to `self.parameters`
        :type span_h: float
        :param span_h: horizontal smoothing length in meters
        :type span_v: float
        :param span_v: vertical smoothing length in meters
        :type output: str
        :param output: file to output stdout to
        :type use_gpu: bool
        :param use_gpu: whether to use GPU acceleration for smoothing. Requires
            GPU compiled binaries and GPU compute node.
        """
        # Redundant to 'base' class but necessary
        if not os.path.exists(input_path):
            unix.mkdir(input_path)

        unix.cd(os.path.join(self.cwd, "DATA"))

        # Copy over only the files that are required. Won't execute if no match
        files = []
        for tag in ["jacobian", "NSPEC_ibool", "x", "y", "z"]:
            files += glob(f"*_{tag}.bin")
        for src in files:
            unix.cp(src=src, dst=input_path)

        super().smooth(input_path=input_path, output_path=output_path,
                       parameters=parameters, span_h=span_h, span_v=span_v,
                       output=output)

    def _initialize_adjoint_traces(self):
        """
        Setup utility: Creates the "adjoint traces" expected by SPECFEM.
        This is only done for the 'base' the Preprocess class.

        .. note::
            Adjoint traces are initialized by writing zeros for all channels.
            Channels actually in use during an inversion or migration will be
            overwritten with nonzero values later on.
        """
        super()._initialize_adjoint_traces()
    
        unix.cd(os.path.join(self.cwd, "traces", "adj"))

        # work around SPECFEM2D's use of different name conventions for
        # regular traces and 'adjoint' traces
        if self.data_format.upper() == "SU":
            files = glob("*SU")
            unix.rename(old="_SU", new="_SU.adj", names=files)
        elif self.data_format.upper() == "ASCII":
            files = glob("*sem?")

            # Get the available extensions, which are named based on unit
            extensions = set([os.path.splitext(_)[-1] for _ in files])
            for extension in extensions:
                unix.rename(old=extension, new=".adj", names=files)

        # SPECFEM2D requires that all components exist even if ununsed
        components = ["x", "y", "z", "p"]

        if self.data_format.upper() == "SU":
            for comp in components:
                src = f"U{self.components[0]}_file_single.su.adj"
                dst = f"U{comp.lower()}s_file_single.su.adj"
                if not os.path.exists(dst):
                    unix.cp(src, dst)
        elif self.data_format.upper() == "ASCII":
            for fid in glob("*.adj"):
                net, sta, cha, ext = fid.split(".")
                for comp in components:
                    # Replace the last value in the channel with new component
                    cha_check = cha[:-1] + comp.upper()
                    fid_check = ".".join([net, sta, cha_check, ext])
                    if not os.path.exists(fid_check):
                        unix.cp(fid, fid_check)
