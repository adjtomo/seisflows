#!/usr/bin/env python3
"""
This class provides utilities for the Seisflows solver interactions with
Specfem2D. It builds upon the base Specfem class which generalizes all solver
interactions with various versions of Specfem.

TODO
    Internal paramater f0 is not currently used. Can we remove or integrate?
"""
import os
from glob import glob

from seisflows.solver.specfem import Specfem
from seisflows.tools import unix
from seisflows.tools.specfem import getpar, setpar


class Specfem2D(Specfem):
    """
    Solver SPECFEM2D
    ----------------
    SPECFEM2D-specific alterations to the base SPECFEM module

    Parameters
    ----------
    :type source_prefix: str
    :param source_prefix: Prefix of source files in path SPECFEM_DATA. Defaults
        to 'SOURCE'
    :type multiples: bool
    :param multiples: set an absorbing top-boundary condition

    Paths
    -----
    ***
    """
    __doc__ = Specfem.__doc__ + __doc__

    def __init__(self, source_prefix="SOURCE", multiples=False, **kwargs):
        """Instantiate a Specfem2D solver interface"""
        super().__init__(source_prefix=source_prefix, **kwargs)

        self.multiples = multiples
        self._f0 = None

        # Define parameters based on material type
        if self.materials.upper() == "ACOUSTIC":
            self._parameters += ["vp"]
        elif self.materials.upper() == "ELASTIC":
            self._parameters += ["vp", "vs"]

    def setup(self):
        """
        Setup the SPECFEM2D solver interface in a SeisFlows workflow
        Append coordinate files to exported model files so that we can use
        them for plotting later
        """
        source_file = os.path.join(self.path.specfem_data, self.source_prefix)
        self._f0 = getpar(key="f0", file=source_file)[1]

        par_file = os.path.join(self.path.specfem_data, "Par_file")
        if self.multiples:
            setpar(key="absorbtop", val=".false.", file=par_file)
        else:
            setpar(key="absorbtop", val=".true.", file=par_file)

        super().setup()

        # Copy in coordinate files to the Model definition so we can plot
        self._export_starting_models(parameters=["x", "z"])

    def smooth(self, input_path, output_path, parameters=None, span_h=None,
               span_v=None, use_gpu=False):
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
            defaults to `self._parameters`
        :type span_h: float
        :param span_h: horizontal smoothing length in meters
        :type span_v: float
        :param span_v: vertical smoothing length in meters
        :type use_gpu: bool
        :param use_gpu: whether to use GPU acceleration for smoothing. Requires
            GPU compiled binaries and GPU compute node.
        """
        unix.cd(os.path.join(self.cwd, self.model_databases))

        # SPECFEM2D requires these files to run the smoother
        files = []
        for tag in ["jacobian", "NSPEC_ibool", "x", "y", "z"]:
            files += glob(f"*_{tag}.bin")
        for src in files:
            unix.cp(src=src, dst=input_path)

        super().smooth(input_path=input_path, output_path=output_path,
                       parameters=parameters, span_h=span_h, span_v=span_v,
                       use_gpu=use_gpu)
