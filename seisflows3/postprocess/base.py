#!/usr/bin/env python
"""
This is the base class for the postprocess functionalities
"""
import os
import sys
import logging

from seisflows3.tools import msg
from seisflows3.tools.wrappers import exists
from seisflows3.config import SeisFlowsPathsParameters

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']


class Base:
    """
    Postprocessing in a Seisflows workflow includes tasks such as
    regularization, smoothing, sharpening, masking and related operations
    on models or gradients
    """
    # Class-specific logger accessed using self.logger
    logger = logging.getLogger(__name__).getChild(__qualname__)

    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings
        """
        pass

    @property
    def required(self):
        """
        A hard definition of paths and parameters required by this class,
        alongside their necessity for the class and their string explanations.
        """
        sf = SeisFlowsPathsParameters()

        # Define the Parameters required by this module
        sf.par("SMOOTH_H", required=False, default=0., par_type=float,
               docstr="Gaussian half-width for horizontal smoothing in units "
                      "of meters. If 0., no smoothing applied")

        sf.par("SMOOTH_V", required=False, default=0., par_type=float,
               docstr="Gaussian half-width for vertical smoothing in units "
                      "of meters")

        sf.par("TASKTIME_SMOOTH", required=False, default=1, par_type=int,
               docstr="Large radii smoothing may take longer than normal "
                      "tasks. Allocate additional smoothing task time "
                      "as a multiple of TASKTIME")

        # Define the Paths required by this module
        sf.path("MASK", required=False, 
                docstr="Directory to mask files for gradient masking")

        return sf

    def check(self, validate=True):
        """
        Checks parameters and paths
        """
        msg.check(type(self))

        if validate:
            self.required.validate()
        
        if PATH.MASK:
            assert exists(PATH.MASK), f"PATH.MASK provided but does not exist"

    def setup(self):
        """
        Placeholder for initialization or setup tasks
        """
        msg.setup(type(self))
        pass

    def write_gradient(self, path):
        """
        Combines contributions from individual sources and material parameters
        to get the gradient, and optionally applies user-supplied scaling

        .. note::
            Because processing operations can be quite expensive, they must be
            run through the HPC system interface; processing does not involve
            embarassingly parallel tasks, we use system.run_single instead of
            system.run

        :type path: str
        :param path: directory from which kernels are read and to which
        gradient is written
        """
        msg.whoami(type(self), prepend="writing gradient from ")

        # Check that the given path exists
        if not exists(path):
            raise FileNotFoundError

        # Run postprocessing as separate job as its computationally intensive
        system.run_single("postprocess", "process_kernels",
                          path=f"{path}/kernels",
                          parameters=solver.parameters,
                          logger=self.logger,
                          scale_tasktime=PAR.TASKTIME_SMOOTH,
                          )

        # Access the gradient information stored in the kernel summation
        gradient = solver.load(f"{path}/kernels/sum", suffix="_kernel")

        # Merge the gradients into a single vector
        gradient = solver.merge(gradient)

        # Convert to absolute perturbations:
        # log dm --> dm (see Eq.13 Tromp et al 2005)
        gradient *= solver.merge(solver.load(f"{path}/model"))

        if PATH.MASK:
            if PAR.VERBOSE:
                print(f"\tMasking gradient")
            # to scale the gradient, users can supply "masks" by exactly
            # mimicking the file format in which models stored
            mask = solver.merge(solver.load(PATH.MASK))

            # While both masking and preconditioning involve scaling the
            # gradient, they are fundamentally different operations:
            # masking is ad hoc, preconditioning is a change of variables;
            # see Modrak & Tromp 2016 GJI
            solver.save(solver.split(gradient), f"{path}/gradient_nomask",
                        parameters=solver.parameters,
                        suffix="_kernel")

            solver.save(solver.split(gradient*mask), f"{path}/gradient",
                        parameters=solver.parameters,
                        suffix="_kernel")
        else:
            solver.save(solver.split(gradient), f"{path}/gradient",
                        parameters=solver.parameters,
                        suffix="_kernel")

    @staticmethod
    def process_kernels(path, parameters, logger):
        """
        Sums kernels from individual sources, with optional smoothing

        Private method because this function needs to be run on system, i.e.,
        called by system.run() or system.run_single(). 

        :type path: str
        :param path: directory containing sensitivity kernels
        :type parameters: list
        :param parameters: solver-defined list of material parameters 
            e.g. ['vp','vs']. This is ideally defined externally by
            PAR.MATERIALS
        :type logger: Logger
        :param logger: Class-specific logging module, log statements pushed
            from this logger will be tagged by its specific module/classname
        """
        # Check if the path exists
        if not exists(path):
            raise FileNotFoundError(f"{path} for postprocess.process_kernels() "
                                    f"not found, exiting.")

        # If specified, smooth the kernels in the vertical and horizontal
        path_sum_nosmooth = os.path.join(path, "sum_nosmooth")
        path_sum = os.path.join(path, "sum")
        if (PAR.SMOOTH_H > 0) or (PAR.SMOOTH_V > 0):
            logger.debug(f"saving unsmoothed, summed kernels to "
                         f"{path_sum_nosmooth}")
            solver.combine(input_path=path, output_path=path_sum_nosmooth,
                           parameters=parameters)

            logger.info(f"smoothing gradient: H={PAR.SMOOTH_H}m, "
                        f"V={PAR.SMOOTH_V}m")
            logger.debug(f"saving smoothed kernels to {path_sum}")
            solver.smooth(input_path=path_sum_nosmooth, output_path=path_sum, 
                          parameters=parameters,
                          span_h=PAR.SMOOTH_H, span_v=PAR.SMOOTH_V)

        # Combine all the input kernels, generating the unscaled gradient
        else:
            logger.debug(f"saving summed kernels to {path_sum}")
            solver.combine(input_path=path, output_path=path_sum,
                           parameters=parameters)
