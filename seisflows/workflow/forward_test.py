"""
Test workflow to see if a new form of seisflows workflow can be used
"""
import os
from seisflows.core import Dict
from seisflows.config import custom_import, config_logger, NAMES
from seisflows.tools.wrappers import load_yaml
from seisflows.tools.specfem import Model



if __name__ == "__main__":
    # Standard SeisFlows Workflow setup block
    # ==========================================================================
    cwd = os.getcwd()
    pars = load_yaml("parameters.yaml")
    logger = config_logger(level=pars.log_level, filename=pars.path_log_file,
                           verbose=pars.verbose)

    # Dynamically create module instances, instantiated with parameters
    classes = [custom_import(name, pars[name.upper()]) for name in NAMES]
    instances = [cls(**pars) for cls in classes]
    # Check that parameters have been set correctly
    for instance in instances:
        instance.check()

    # Distribute instances to their respective namesakes
    system, preprocess, solver, postprocess, optimize, workflow = instances
    # ==========================================================================


    # Begin workflow
    logger.info("Starting forward simulation workflow")


    m = Model(pars.path_model_init)
    m.write(path=os.path.join(paths.eval_func, "model"))

    system.run(solver.eval_func)

