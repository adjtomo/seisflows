#!/usr/bin/env python3
"""
This class provides the core utilities for the SeisFlows postprocessing
functionalities, including kernel/gradient smoothing and masking as well as
kernel summation
"""
import os

from seisflows import logger
from seisflows.tools.specfem import Model


class Default:
    """
    Postprocessing in a Seisflows workflow includes tasks such as
    regularization, smoothing, sharpening, masking and related operations
    on models or gradients
    """
    def __init__(self, smooth_h=0., smooth_v=0., tasktime_smooth=1,
                 path_postprocess=None, path_mask=None, **kwargs):
        """
        Establish Postprocessing parameters

        :type smooth_h: float
        :param smooth_h: Gaussian half-width for horizontal smoothing in units
            of meters. If 0., no smoothing applied
        :type smooth_h: float
        :param smooth_v: Gaussian half-width for vertical smoothing in units
            of meters.
        :type tasktime_smooth: float
        :param tasktime_smooth: Large radii smoothing may take longer than
            normal tasks. Allocate additional smoothing task time as a multiple
            of system.tasktime
        :type path_postprocess: str
        :param path_postprocess: scratch path to perform all postprocessing
            tasks such as smoothing, kernel combination etc.
        :type path_mask: str
        :param path_mask: Directory to mask files for gradient masking. Format
            of the mask files MUST match the format of the input model.
        """
        super().__init__()

        self.smooth_h = smooth_h
        self.smooth_v = smooth_v
        self.tasktime_smooth = tasktime_smooth
        self.path = path_postprocess or \
                    os.path.join(os.getcwd(), "scratch", "evalgrad")
        self.path_mask = path_mask

    def check(self):
        """
        Checks parameters and paths
        """
        if self.path_mask:
            assert os.path.exists(self.path_mask), \
                "`postprocess.path_mask` provided but does not exist"

    def setup(self):
        """
        Setup tasks
        """
        pass

    def finalize(self):
        """
        Finalization tasks
        """
        pass

    def scale_gradient(self):
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
        # Access the gradient information stored in as kernel files
        model = Model(path=os.path.join(self.path, "model"))
        gradient = Model(path=os.path.join(self.path, "kernels", "sum"))

        # Merge to vector and convert to absolute perturbations:
        # log dm --> dm (see Eq.13 Tromp et al 2005)
        gradient.vector *= model.vector

        if self.path_mask:
            logger.info(f"masking gradient")
            # to scale the gradient, users can supply "masks" by exactly
            # mimicking the file format in which models are stored
            mask = Model(self.path_mask)

            gradient.write(path=os.path.join(self.path, "gradient_nomask"))

            gradient.vector *= mask.vector

        return gradient

    def sum_smooth_kernels(self, solver):
        """
        Sums kernels from individual sources, with optional smoothing

        .. note::
            This function needs to be run on system, i.e., called by
            system.run(single=True)

        :type solver: solver instance
        :param solver: SeisFlows solver which will be used for its combine and
            smooth functions
        """
        # If specified, smooth the kernels in the vertical and horizontal and
        # save both (summed, summed+smoothed) to separate output directories
        path_kernel = os.path.join(self.path, "kernels")
        path_sum_nosmooth = os.path.join(path_kernel, "sum_nosmooth")
        path_sum = os.path.join(path_kernel, "sum")

        if (self.smooth_h > 0) or (self.smooth_v > 0):
            logger.debug(f"saving un-smoothed and summed kernels to:\n"
                         f"{path_sum_nosmooth}")
            solver.combine(input_path=path_kernel, output_path=path_sum_nosmooth)

            logger.info(f"smoothing gradient: "
                        f"H={self.smooth_h}m; V={self.smooth_v}m")
            logger.debug(f"saving smoothed kernels to:\n{path_sum}")
            solver.smooth(input_path=path_sum_nosmooth, output_path=path_sum, 
                          span_h=self.smooth_h, span_v=self.smooth_v)

        # Combine all the input kernels, generating the unscaled gradient
        else:
            logger.debug(f"saving summed kernels to:\n{path_sum}")
            solver.combine(input_path=path_kernel, output_path=path_sum)
