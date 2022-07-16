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
from seisflows.workflow.forward import Forward
from seisflows.tools import msg, unix
from seisflows.tools.specfem import Model


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
        """
        Init is used to instantiate global parameters defined by the input
        parameter file.
        """
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
                self.generate_misfit_kernels,
                self.postprocess_kernels
                ]

    def setup(self):
        """
        Override the Forward.setup() method to include the postprocessing
        module used for kernel/gradient manipulation
        """
        super().setup()
        self.postprocess = self._modules.postprocess

    def generate_misfit_kernels(self):
        """System wrapper for running adjoint simulations"""
        logger.msg.mnr("GENERATING MISFIT KERNELS")
        self.system.run(self.generate_misfit_kernel)

    def generate_misfit_kernel(self):
        """
        Performs adjoint simulations for a single given event. File manipulation
        to ensure kernels are discoverable by other modules
        """
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

    def postprocess_kernels(self):
        """System wrapper for postprocess kernels. Run with single"""
        self.system.run(self._postprocess_kernels, single=True)

    def _postprocess_kernels(self):
        """
        System-run wrapper for postprocess.process_kernels which is meant to
        sum and smooth all individual event kernels
        """
        # Combine kernels into a single volumentric quantity
        self.solver.combine(
            input_path=os.path.join(self.path.eval_grad, "kernels"),
            output_path=os.path.join(self.path.eval_grad, "sum")
        )

        if self.solver.smooth_h > 0. or self.solver.smooth_v > 0.:
            # Make a distinction that we have a pre- and post-smoothed sum
            unix.mv(src=os.path.join(self.path.eval_grad, "sum_nosmooth"))

            self.solver.smooth(
                input_path=os.path.join(self.path.eval_grad, "sum_nosmooth"),
                output_path=os.path.join(self.path.eval_grad, "sum")
            )
