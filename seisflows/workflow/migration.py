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
import numpy as np
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

    def __init__(self, path_mask=None, export_gradient=True,
                 export_kernels=False, **kwargs):
        """
        Instantiate Migration-specific parameters
        """
        super().__init__(**kwargs)

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
        return [self.generate_synthetic_data,
                self.evaluate_initial_misfit,
                self.run_adjoint_simulations,
                self.postprocess_event_kernels,
                self.evaluate_gradient_from_kernels,
                self.finalize_iteration,
                ]

    def run_adjoint_simulations(self, **kwargs):
        """
        Performs adjoint simulations for all events. Some additional file 
        naming to ensure kernels are discoverable by other modules. 

        .. note::

            This is a relatively simple function, but it keeps the approach 
            general by allowing other workflows to include pre- and post-
            processing tasks, or to overwrite the adjoint simulation task
        """
        self.system.run([self._run_adjoint_simulation_single], **kwargs)

    def _run_adjoint_simulation_single(self, save_kernels=None, 
                                       export_kernels=None, **kwargs):
        """
        Run an adjoint simulation for a single source. Allow saving kernels by 
        moving them out of the run directory to another location. Allow 
        exporting kernels by copying them to the output directory.

        .. note::

            Must be run by system.run() so that solvers are assigned
            individual task ids/working directories.

        .. note:: 

            see solver.specfem.adjoint_simulation() for full 
            detailed list of input parameters

        :type save_kernels: str
        :param save_kernels: path to a directory where kernels created by the 
            adjoint simulation are moved to for further use in the workflow
            Defaults to saving kernels in `scratch/eval_grad/kernels/<source>`
        :type export kernels: str
        :param export_kernels: path to a directory where kernels are copied for
            more permanent storage, where they will not be wiped by `clean` or
            `restart`. User parameter `export_kernels` must be set `True`.
        """
        # Set default value for `export_kernels` or take program default
        if self.export_kernels:
            if export_kernels is None:
                export_kernels = os.path.join(self.path.output, "kernels",
                                              self.solver.source_name)
        else:
            export_kernels = False

        # Set default value for `save_kernels` or take programmed default
        if save_kernels is None:
            save_kernels = os.path.join(self.path.eval_grad, "kernels",
                                        self.solver.source_name, "")

        logger.info(f"running adjoint simulation for source "
                    f"{self.solver.source_name}")

        # Run adjoint simulations on system. Make kernels discoverable in
        # path `eval_grad`. Optionally export those kernels
        self.solver.adjoint_simulation(save_kernels=save_kernels,
                                       export_kernels=export_kernels,
                                       **kwargs)

    def postprocess_event_kernels(self):
        """
        Combine/sum NTASK event kernels into a single volumetric kernel and
        then (optionally) smooth the output misfit kernel by convolving with
        a 3D Gaussian function with user-defined horizontal and vertical
        half-widths, or by using the laplacian smoothing operator.

        .. note::

            If you hit a floating point error during the smooth operation, 
            your kernels may be zero due to something going awry in the
            misfit quantification or adjoint simulations.
        """
        def mask_source_event_kernels(**kwargs):
            """
            Mask source region by combining binary source mask with kernel.
            This feature is only available in SPECFEM3D_GLOBE and only turned on
            if solver.mask_source is set to True, otherwise it will be skipped.

            This uses the Model class because SPECFEM does not have an internal
            function for multiplying files (only for adding/subtracting)
            """
            # Only trigger this function if the Solver saved source mask files
            mask_path = os.path.join(self.path.eval_grad, "mask_source", 
                                     self.solver.source_names[0])
            if not glob(os.path.join(mask_path, "*")):
                logger.debug("no source mask files found, skipping source mask")
                return
            
            logger.info("masking source region in event kernels")
            for src in self.solver.source_names:
                logger.debug(f"mask source {src}")
                path = os.path.join(self.path.eval_grad, "mask_source", src)

                # Mask vector [0, 1], where values <1 are near source
                mask_model = Model(path=path, parameters=["mask_source"],
                                   regions=self.solver._regions)
                
                # Gaussian mask sometimes does not sufficiently suppress source 
                # region so we modify the amplitude 
                maskv = mask_model.vector
                if self.solver.scale_mask_region:
                    logger.info(f"scaling source mask by "
                                f"{self.solver.scale_mask_region}")
                    maskv[maskv < 1] *= self.solver.scale_mask_region  

                # Need to expand vector the length of the event kernel vector
                maskv = np.tile(maskv, len(self.solver._parameters))
                
                # Now we apply the mask to the event kernel which contains all
                # parameters we are updating in our inversion
                event_kernel = Model(
                    path=os.path.join(self.path.eval_grad, "kernels", src), 
                    parameters=[f"{par}_kernel" for par in 
                                self.solver._parameters],
                    regions=self.solver._regions
                    )
                event_kernel.update(vector=event_kernel.vector * maskv)

                # Overwrite existing event kernel files with the masked version
                event_kernel.write(
                    path=os.path.join(self.path.eval_grad, "kernels", src)
                    )

        def combine_event_kernels(**kwargs):
            """
            Combine individual event kernels into a misfit kernel for each
            parameter defined by the solver
            """
            # Input paths are the kernels generated by each of the sources
            input_paths = [os.path.join(self.path.eval_grad, "kernels", src) for
                           src in self.solver.source_names]

            logger.info("combining event kernels into single misfit kernel")

            # Parameters to combine are the kernels, which follow the 
            # naming convention {par}_kernel
            parameters = [f"{par}_kernel" for par in self.solver._parameters]

            self.solver.combine(
                input_paths=input_paths,
                output_path=os.path.join(self.path.eval_grad, "misfit_kernel"),
                parameters=parameters
            )

        def smooth_misfit_kernel(**kwargs):
            """
            Smooth the misfit kernel using the underlying Solver smooth function
            """
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

        self.system.run([mask_source_event_kernels, combine_event_kernels, 
                         smooth_misfit_kernel], single=True)

    def evaluate_gradient_from_kernels(self):
        """
        Generates the 'gradient' from the 'misfit kernel'. This involves
        scaling the gradient by the model vector (log dm --> dm) and applying
        an optional mask function to the gradient.

        :raises SystemError: if the gradient vector is zero, which means that
            none of the kernels returned usable values.
        """
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
        logger.info("scaling gradient to absolute model perturbations")
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
            dst = os.path.join(self.solver.path.output, "gradient")
            unix.cp(src, dst)

        # Last minute check to see if the gradient is 0. If so then the adjoint
        # simulations returned no kernels and that probably means something 
        # went wrong
        if sum(gradient.vector) == 0:
            logger.critical(msg.cli(
                "Gradient vector 'g' is 0, meaning none of the kernels from "
                "the adjoint simulations returned. Please check your gradient "
                "and waveform misfits. ", border="=",
                header="optimization gradient error")
            )
            sys.exit(-1)

