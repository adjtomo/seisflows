"""
Test workflow to see if a new form of seisflows workflow can be used
"""
import os
from glob import glob
from seisflows import logger
from seisflows.config import import_seisflows
from seisflows.tools import msg

# Standard SeisFlows setup, makes modules global variables to the workflow
pars, modules = import_seisflows()
system, preprocess, solver, postprocess, optimize = modules


def evaluate_objective_function(path_model):
    """
    Performs forward simulation for a single given event. Also evaluates the
    objective function and writes residuals and adjoint sources for later tasks.

    .. note::
        if PAR.PREPROCESS == None, will not perform misfit quantification

    .. note::
        Must be run by system.run() so that solvers are assigned individual
        task ids/ working directories.
    """
    # Run the forward simulation with the given input model
    solver.import_model(path_model=path_model)
    solver.forward_simulation(
        save_traces=os.path.join(solver.cwd, "traces", "syn"),
        export_traces=os.path.join(solver.path.output,
                                   solver.source_name, "syn")
    )

    # Perform data-synthetic misfit quantification
    if preprocess is not None:
        preprocess.quantify_misfit(
            observed=solver.data_filenames(choice="obs"),
            synthetics=solver.data_filenames(choice="syn"),
            output=os.path.join(solver.cwd, "traces", "adj")
        )


if __name__ == "__main__":
    # Begin the forward simulation workflow
    logger.info(msg.mjr("Starting forward simulation workflow"))

    for module in modules:
        module.setup()

    # Run objective function evaluation NTASK times
    logger.info(msg.sub("EVALUATING OBJECTIVE FUNCTION"))
    system.run(evaluate_objective_function, path_model=pars.path_model_init)

    logger.info(msg.mjr("Finished forward simulation workflow"))

