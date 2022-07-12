"""
Test workflow to see if a new form of seisflows workflow can be used
"""
import os
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
    """
    if system.taskid == 0:
        logger.info(msg.sub("EVALUATING OBJECTIVE FUNCTION"))

    # Run the forward simulation with the given input model
    solver.import_model(path=path_model)
    solver.forward_simulation(
        output_seismograms=os.path.join(solver.cwd, "traces", "syn")
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
    system.run(evaluate_objective_function, path_model=pars.path_model_init,
               suffix="new", system=system)

    logger.info(msg.mjr("Finished forward simulation workflow"))

