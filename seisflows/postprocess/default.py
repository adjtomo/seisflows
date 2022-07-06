#!/usr/bin/env python3
"""
This class provides the core utilities for the SeisFlows postprocessing
functionalities, including kernel/gradient smoothing and masking as well as
kernel summation
"""
import os
import sys

from seisflows.core import Base
from seisflows.tools import msg


class Default(Base):
    """
    Postprocessing in a Seisflows workflow includes tasks such as
    regularization, smoothing, sharpening, masking and related operations
    on models or gradients
    """
    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings
        """
        super().__init__()

        self.required.par(
            "SMOOTH_H", required=False, default=0., par_type=float,
            docstr="Gaussian half-width for horizontal smoothing in units of "
                   "meters. If 0., no smoothing applied"
        )
        self.required.par(
            "SMOOTH_V", required=False, default=0., par_type=float,
            docstr="Gaussian half-width for vertical smoothing in units of "
                   "meters"
        )
        self.required.par(
            "TASKTIME_SMOOTH", required=False, default=1, par_type=int,
            docstr="Large radii smoothing may take longer than normal tasks. "
                   "Allocate additional smoothing task time as a multiple of "
                   "TASKTIME"
        )
        # Define the Paths required by this module
        self.required.path(
            "MASK", required=False, docstr="Directory to mask files for "
                                           "gradient masking"
        )

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        super().check(validate=validate)
        
        if self.path.MASK:
            assert os.path.exists(self.path.MASK), \
                f"PATH.MASK provided but does not exist"

    def setup(self):
        """
        Setup tasks
        """
        super().setup()

    def finalize(self):
        """
        Finalization tasks
        """
        super().finalize()

    def scale_gradient(self, input_path):
        """
        Combines contributions from individual sources and material parameters
        to get the gradient, and optionally applies user-supplied scaling

        .. note::
            Because processing operations can be quite expensive, they must be
            run through the HPC system interface; processing does not involve
            embarassingly parallel tasks, we use run(single=True)

        :type input_path: str
        :param input_path: directory from which kernels are read and to which
            gradient is written. Should probably point to PATH.GRAD
        :rtype: np.array
        :return: scaled gradient as a vector
        """
        solver = self.module("solver")

        # Postprocess file structure defined here once-and-for-all
        path_grad_nomask = os.path.join(input_path, "gradient_nomask")
        path_model = os.path.join(input_path, "model")
        path_kernels_sum = os.path.join(input_path, "kernels", "sum")

        # Access the gradient information stored in as kernel files
        gradient = solver.load(path_kernels_sum, suffix="_kernel")

        # Merge to vector and convert to absolute perturbations:
        # log dm --> dm (see Eq.13 Tromp et al 2005)
        gradient = solver.merge(gradient)
        gradient *= solver.merge(solver.load(path_model))

        if self.path.MASK:
            self.logger.info(f"masking gradient")
            # to scale the gradient, users can supply "masks" by exactly
            # mimicking the file format in which models are stored
            mask = solver.merge(solver.load(self.path.MASK))

            # While both masking and preconditioning involve scaling the
            # gradient, they are fundamentally different operations:
            # masking is ad hoc, preconditioning is a change of variables;
            # For more info, see Modrak & Tromp 2016 GJI
            solver.save(solver.split(gradient), path=path_grad_nomask,
                        suffix="_kernel")
            gradient *= mask

        return gradient

    def sum_smooth_kernels(self, path_grad):
        """
        Sums kernels from individual sources, with optional smoothing

        .. note::
            This function needs to be run on system, i.e., called by
            system.run(single=True)

        :type path_grad: str
        :param path_grad: directory containing sensitivity kernels in the
            scratch directory to be summed and smoothed. Output summed and
            summed + smoothed kernels will be saved here as well.
        """
        solver = self.module("solver")

        # If specified, smooth the kernels in the vertical and horizontal and
        # save both (summed, summed+smoothed) to separate output directories
        kernel_path = os.path.join(path_grad, "kernels")

        path_sum_nosmooth = os.path.join(kernel_path, "sum_nosmooth")
        path_sum = os.path.join(kernel_path, "sum")

        if (self.par.SMOOTH_H > 0) or (self.par.SMOOTH_V > 0):
            self.logger.debug(f"saving un-smoothed and summed kernels to:\n"
                              f"{path_sum_nosmooth}")
            solver.combine(input_path=kernel_path,
                           output_path=path_sum_nosmooth)

            self.logger.info(f"smoothing gradient: H={self.par.SMOOTH_H}m, "
                             f"V={self.par.SMOOTH_V}m")
            self.logger.debug(f"saving smoothed kernels to:\n{path_sum}")
            solver.smooth(input_path=path_sum_nosmooth, output_path=path_sum, 
                          span_h=self.par.SMOOTH_H, span_v=self.par.SMOOTH_V)

        # Combine all the input kernels, generating the unscaled gradient
        else:
            self.logger.debug(f"saving summed kernels to:\n{path_sum}")
            solver.combine(input_path=kernel_path, output_path=path_sum)
