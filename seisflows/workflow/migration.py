#!/usr/bin/env python3
"""
Seismic migration performs a 'time-reverse migration', or backprojection.
In the terminology of seismic imaging, we are running a forward and adjoint
simulation to derive the gradient of the objective function. This workflow
sets up the machinery to derive a scaled, smoothed gradient from an initial
model
"""
import os

from seisflows import logger
from seisflows.config import import_seisflows
from seisflows.workflow.forward import Forward
from seisflows.tools import msg, unix


class Migration(Forward):
    """
    [workflow.migration] Run forward and adjoint solver to produce
    event-dependent misfit kernels. Sum and postprocess kernels to produce
    gradient. In seismic exploration this is 'reverse time migration'.

    .. warning::
        Misfit kernels require large amounts of disk space for storage.
        Setting `export_kernel`==True when PAR.NTASK is large and model files
        are large may lead to large file overhead.

    :type export_gradient: bool
    :param export_gradient: export the gradient after it has been generated
        in the scratch directory. If False, gradient can be discarded from
        scratch at any time in the workflow
    :type export_kernels: bool
    :param export_kernels: export each sources event kernels after they have
        been generated in the scratch directory. If False, gradient can be
        discarded from scratch at any time in the workflow
    """
    __doc__ = Forward.__doc__ + __doc__

    def __init__(self, _modules=None, export_gradient=False,
                 export_kernels=False, **kwargs):
        """Instantiate Migration-specific parameters"""
        super().__init__(**kwargs)

        self._modules = _modules
        self.export_gradient = export_gradient
        self.export_kernels = export_kernels

        # Overwriting base class required modules list
        self._required_modules = ["system", "solver", "preprocess",
                                  "postprocess"]

        # Empty module variables that should be filled in by setup
        self.postprocess = None

    @property
    def task_list(self):
        """
        USER-DEFINED TASK LIST. This property defines a list of class methods
        that take NO INPUT and have NO RETURN STATEMENTS. This defines your
        linear workflow, i.e., these tasks are to be run in order from start to
        finish to complete a workflow.

        This excludes 'check' (which is run during 'import_seisflows') and
        'setup' which should be run separately

        .. note::
            For workflows that require an iterative approach (e.g. inversion),
            this task list will be looped over, so ensure that any setup and
            teardown tasks (run once per workflow, not once per iteration) are
            not included.

        :rtype: list
        :return: list of methods to call in order during a workflow
        """
        return [self.evaluate_initial_misfit,
                self.run_adjoint_simulations,
                self.generate_misfit_kernel
                ]

    def setup(self):
        """
        Override the Forward.setup() method to include the postprocessing
        module used for kernel/gradient manipulation
        """
        super().setup()
        self.postprocess = self._modules.postprocess

    def run_adjoint_simulations(self):
        """
        Performs adjoint simulations for a single given event. File manipulation
        to ensure kernels are discoverable by other modules
        """
        def run_adjoint_simulation():
            """Adjoint simulation function to be run by system.run()"""
            if self.export_kernels:
                export_kernels = os.path.join(self.path.output, "kernels",
                                              self.solver.source_name)
            else:
                export_kernels = False

            # Run adjoint simulations on system. Make kernels discoverable in
            # path `eval_grad`. Optionally export those kernels
            self.solver.adjoint_simulation(
                save_kernels=os.path.join(self.path.eval_grad, "kernels",
                                          self.solver.source_name),
                export_kernels=export_kernels
            )

        logger.msg.mnr("running adjoint simulations to generate event kernels")
        self.system.run(run_adjoint_simulation)

    def generate_misfit_kernel(self):
        """
        Combine/sum NTASK event kernels into a single volumetric kernel and
        then (optionally) smooth the output misfit kernel by convolving with
        a 3D Gaussian function with user-defined horizontal and vertical
        half-widths.
        """
        def combine_event_kernels():
            """Combine event kernels into a misfit kernel"""
            self.solver.combine(
                input_path=os.path.join(self.path.eval_grad, "kernels"),
                output_path=os.path.join(self.path.eval_grad, "sum")
            )

        def smooth_misfit_kernel():
            """Smooth the output misfit kernel """
            if self.solver.smooth_h > 0. or self.solver.smooth_v > 0.:
                # Make a distinction that we have a pre- and post-smoothed sum
                unix.mv(src=os.path.join(self.path.eval_grad, "sum_nosmooth"))

                self.solver.smooth(
                    input_path=os.path.join(self.path.eval_grad, "sum_nosmooth"),
                    output_path=os.path.join(self.path.eval_grad, "sum")
                )

        logger.msg.mnr("postprocessing (summing/smoothing) event kernels")
        self.system.run([combine_event_kernels, smooth_misfit_kernel],
                        single=True)


if __name__ == "__main__":
    # Standard SeisFlows setup, makes modules global variables to the workflow
    pars, modules = import_seisflows()

    logger.info(msg.mjr("Starting migration workflow"))

    workflow = Migration(modules, **pars)
    workflow.check()
    workflow.setup()
    workflow.run()

    logger.info(msg.mjr("Finished migration workflow"))