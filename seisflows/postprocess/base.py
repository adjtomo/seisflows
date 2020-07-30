#!/usr/bin/env python
"""
This is the base class for the postprocess functionalities
"""
import sys

from seisflows.tools.tools import exists

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
    def __init__(self):
        """
        These parameters should not be set by __init__!
        Attributes are just initialized as NoneTypes for clarity and docstrings
        """

    @staticmethod
    def check():
        """
        Checks parameters and paths
        """
        if "SMOOTH_H" not in PAR:
            setattr(PAR, "SMOOTH_H", 0.)
        
        if "SMOOTH_V" not in PAR:
            setattr(PAR, "SMOOTH_V", 0.)

        if "TASKTIME_SMOOTH" not in PAR:
            setattr(PAR, "TASKTIME_SMOOTH", 1)

        if "MASK" not in PATH:
            setattr(PATH, "MASK", None)

        if PATH.MASK:
            assert exists(PATH.MASK)

    def setup(self):
        """
        Placeholder for initialization or setup tasks
        """
        pass

    @staticmethod
    def process_kernels(path, parameters):
        """
        Sums kernels from individual sources, with optional smoothing
        This should be run in paralell

        :type path: str
        :param path: directory containing sensitivity kernels
        :type parameters: list
        :param parameters: material parameters e.g. ['vp','vs']
        """
        # Check if the path exists
        if not exists(path):
            raise FileNotFoundError

        # If specified, smooth the kernels in the vertical and horizontal
        if PAR.SMOOTH_H > 0:
            solver.combine(input_path=path, output_path=f"{path}/sum_nosmooth",
                           parameters=parameters)
            print(f"\tSmoothing gradient: H={PAR.SMOOTH_H}, V={PAR.SMOOTH_V}")
            solver.smooth(input_path=f"{path}/sum_nosmooth",
                          output_path=f"{path}/sum", parameters=parameters,
                          span_h=PAR.SMOOTH_H, span_v=PAR.SMOOTH_V)
        else:
            solver.combine(input_path=path, output_path=f"{path}/sum",
                           parameters=parameters)

    @staticmethod
    def write_gradient(path):
        """
        Combines contributions from individual sources and material parameters
        to get the gradient, and optionally applies user-supplied scaling

        Note:
            Because processing operations can be quite expensive, they must be
            run through the HPC system interface; processing does not involve
            embarassingly parallel tasks, we use system.run_single instead of
            system.run

        :type path: str
        :param path: directory from which kernels are read and to which
        gradient is written
        """
        # Check that the given path exists
        if not exists(path):
            raise FileNotFoundError

        # Run postprocessing on the cluster
        system.run_single("postprocess", "process_kernels",
                          path=f"{path}/kernels",
                          scale_tasktime=PAR.TASKTIME_SMOOTH,
                          parameters=solver.parameters)

        # Access the gradient information stored in the kernel summation
        gradient = solver.load(f"{path}/kernels/sum", suffix="_kernel")

        # Merge the gradients into a single vector
        gradient = solver.merge(gradient)

        # Convert to absolute perturbations:
        # log dm --> dm (see Eq.13 Tromp et al 2005)
        gradient *= solver.merge(solver.load(f"{path}/model"))

        if PATH.MASK:
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

