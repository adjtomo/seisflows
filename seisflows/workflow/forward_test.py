"""
Test workflow to see if a new form of seisflows workflow can be used
"""
import os
from glob import glob
from seisflows import logger
from seisflows.config import import_seisflows
from seisflows.tools import msg
from seisflows.tools.utils import log_status
from seisflows.config import config_logger

# Standard SeisFlows setup, makes modules global variables to the workflow
pars, modules = import_seisflows()
# config_logger(level=pars.log_level, filename=pars.path_log_
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
    if system.taskid == 0:
        logger.info(msg.sub("EVALUATING OBJECTIVE FUNCTION"))

    # Run the forward simulation with the given input model
    solver.import_model(path_model=path_model)
    solver.forward_simulation(
        save_traces=os.path.join(solver.cwd, "traces", "syn"),
        export_traces=os.path.join(pars.path_output, solver.source_name, "syn")
    )

    # Perform data-synthetic misfit quantification
    if preprocess is not None:
        preprocess.quantify_misfit(
            observed=solver.data_filenames(choice="obs"),
            synthetics=solver.data_filenames(choice="syn"),
            output=os.path.join(solver.cwd, "traces", "adj")
        )


@log_status
def run_forward_simulation(path_model):
    """

    """
    # Run the forward simulation with the given input model
    solver.import_model(path_model=path_model)
    solver.forward_simulation(
        save_traces=os.path.join(solver.cwd, "traces", "syn"),
        export_traces=os.path.join(pars.path_output, solver.source_name, "syn")
    )


@log_status
def quantify_misfit():
    """
    Quantify the data-synthetic misfit, write residuals to scratch for use in
    potential optimization, generate adjoint sources required for adjoint
    simulations
    """
    preprocess.quantify_misfit(
        observed=solver.data_filenames(choice="obs"),
        synthetics=solver.data_filenames(choice="syn"),
        write_adjsrcs=os.path.join(solver.cwd, "traces", "adj"),
        write_residuals=os.path.join(pars.path_eval_grad, solver.source_name)
    )
    preprocess.sum_residuals(files=glob(os.path.join(pars.path_eval_grad, "*")))


if __name__ == "__main__":
    # Begin the forward simulation workflow
    logger.info(msg.mjr("Starting forward simulation workflow"))

    for module in modules:
        module.setup()

    # Run objective function evaluation NTASK times
    system.run(evaluate_objective_function, path_model=pars.path_model_init,
               suffix="new", system=system)

    logger.info(msg.mjr("Finished forward simulation workflow"))

