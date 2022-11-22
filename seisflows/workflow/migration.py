#!/usr/bin/env python3
"""
Seismic migration performs a 'time-reverse migration', or backprojection.
In the terminology of seismic imaging, we are running a forward and adjoint
simulation to derive the gradient of the objective function. This workflow
sets up the machinery to derive a scaled, smoothed gradient from an initial
model

.. warning::
    Misfit kernels require large amounts of disk space for storage.
    Setting `export_kernel`==True when PAR.NTASK is large and model files
    are large may lead to large file overhead.

.. note::
    Migration workflow includes an option to mask the gradient. While both
    masking and preconditioning involve scaling the gradient, they are
    fundamentally different operations: masking is ad hoc, preconditioning
    is a change of variables; For more info, see Modrak & Tromp 2016 GJI
"""
import os
import sys
import shutil
from glob import glob
from seisflows import logger
from seisflows.tools import msg, unix
from seisflows.tools.model import Model
from seisflows.workflow.forward import Forward


class Migration(Forward):
    """
    Migration Workflow
    ------------------
    Run forward and adjoint solver to produce event-dependent misfit kernels.
    Sum and postprocess kernels to produce gradient. In seismic exploration
    this is 'reverse time migration'.

    Parameters
    ----------
    :type export_gradient: bool
    :param export_gradient: export the gradient after it has been generated
        in the scratch directory. If False, gradient can be discarded from
        scratch at any time in the workflow
    :type export_kernels: bool
    :param export_kernels: export each sources event kernels after they have
        been generated in the scratch directory. If False, gradient can be
        discarded from scratch at any time in the workflow

    Paths
    -----
    :type path_mask: str
    :param path_mask: optional path to a masking function which is used to
        mask out or scale parts of the gradient. The user-defined mask must
        match the file format of the input model (e.g., .bin files).
    ***
    """
    __doc__ = Forward.__doc__ + __doc__

    def __init__(self, modules=None, path_mask=None, export_gradient=True,
                 export_kernels=False, **kwargs):
        """
        Instantiate Migration-specific parameters

        :type modules: list
        :param modules: list of sub-modules that will be established as class
            attributes by the setup() function. Should not need to be set by the
            user
        """
        super().__init__(**kwargs)

        self._modules = modules
        self.export_gradient = export_gradient
        self.export_kernels = export_kernels

        self.path["mask"] = path_mask

        # Overwriting base class required modules list
        self._required_modules = ["system", "solver", "preprocess"]

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
                self.postprocess_event_kernels,
                self.evaluate_gradient_from_kernels
                ]

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

            logger.info(f"running adjoint simulation for source "
                        f"{self.solver.source_name}")
            # Run adjoint simulations on system. Make kernels discoverable in
            # path `eval_grad`. Optionally export those kernels
            self.solver.adjoint_simulation(
                save_kernels=os.path.join(self.path.eval_grad, "kernels",
                                          self.solver.source_name, ""),
                export_kernels=export_kernels
            )

        logger.info(msg.mnr("EVALUATING EVENT KERNELS W/ ADJOINT SIMULATIONS"))
        self.system.run([run_adjoint_simulation])

    def postprocess_event_kernels(self):
        """
        Combine/sum NTASK event kernels into a single volumetric kernel and
        then (optionally) smooth the output misfit kernel by convolving with
        a 3D Gaussian function with user-defined horizontal and vertical
        half-widths.
        """
        def combine_event_kernels():
            """Combine event kernels into a misfit kernel"""
            logger.info("combining event kernels into single misfit kernel")
            self.solver.combine(
                input_path=os.path.join(self.path.eval_grad, "kernels"),
                output_path=os.path.join(self.path.eval_grad, "misfit_kernel")
            )

        def smooth_misfit_kernel():
            """Smooth the output misfit kernel"""
            if self.solver.smooth_h > 0. or self.solver.smooth_v > 0.:
                logger.info(
                    f"smoothing misfit kernel: "
                    f"H={self.solver.smooth_h}; V={self.solver.smooth_v}"
                )
                # Make a distinction that we have a pre- and post-smoothed kern.
                unix.mv(
                    src=os.path.join(self.path.eval_grad, "misfit_kernel"),
                    dst=os.path.join(self.path.eval_grad, "mk_nosmooth")
                )
                self.solver.smooth(
                    input_path=os.path.join(self.path.eval_grad, "mk_nosmooth"),
                    output_path=os.path.join(self.path.eval_grad,
                                             "misfit_kernel")
                )

        # Make sure were in a clean scratch eval_grad directory
        tags = ["misfit_kernel", "mk_nosmooth"]
        for tag in tags:
            scratch_path = os.path.join(self.path.eval_grad, tag)
            if os.path.exists(scratch_path):
                shutil.rmtree(scratch_path)

        logger.info(msg.mnr("GENERATING/PROCESSING MISFIT KERNEL"))
        self.system.run([combine_event_kernels, smooth_misfit_kernel],
                        single=True)

    def evaluate_gradient_from_kernels(self):
        """
        Generates the 'gradient' from the 'misfit kernel'. This involves
        scaling the gradient by the model vector (log dm --> dm) and applying
        an optional mask function to the gradient.
        """
        logger.info("scaling gradient to absolute model perturbations")

        # Check that kernel files exist before attempting to manipulate
        misfit_kernel_path = os.path.join(self.path.eval_grad, "misfit_kernel")
        if not glob(os.path.join(misfit_kernel_path, "*")):
            logger.critical(msg.cli(
                "directory 'scratch/eval_grad/misfit_kernel' is empty but "
                "should contain summed kernels. Please check "
                "'scratch/solver/mainsolver' log files to see if the "
                "`xcombine` and `xsmooth` operations completed successfully", 
                header="missing kernels error", border="=")
                )
            sys.exit(-1)
        # Read from files, only pick up regions defined by solver
        gradient = Model(path=misfit_kernel_path, regions=self.solver._regions)

        # Set model: we only need to access parameters which will be updated
        # Assuming that the model in the solver also generated the kernels
        dst = os.path.join(self.path.eval_grad, "model")
        unix.rm(dst)
        unix.mkdir(dst)
        for src in self.solver.model_files:
            unix.ln(src, dst=os.path.join(dst, os.path.basename(src)))

        # Read in new model files that will have been generated by `setup` 
        # or by optimization library
        model = Model(path=dst, parameters=self.solver._parameters,
                      regions=self.solver._regions)

        # Merge to vector and convert to absolute perturbations:
        # log dm --> dm (see Eq.13 Tromp et al 2005)
        gradient.update(vector=gradient.vector * model.vector)
        gradient.write(path=os.path.join(self.path.eval_grad, "gradient"))

        # Apply an optional mask to the gradient
        if self.path.mask:
            logger.info("applying mask function to gradient")
            mask = Model(path=self.path.mask)
            unix.mv(src=os.path.join(self.path.eval_grad, "gradient"),
                    dst=os.path.join(self.path.eval_grad, "gradient_nomask"))

            gradient.update(vector=gradient.vector * mask.vector)
            gradient.write(path=os.path.join(self.path.eval_grad, "gradient"))

        # Export gradient to disk
        if self.export_gradient:
            logger.info("exporting gradient to disk")
            src = os.path.join(self.path.eval_grad, "gradient")
            dst = os.path.join(self.path.output, "gradient")
            unix.cp(src, dst)

