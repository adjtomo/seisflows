#!/bin/env python

import sys

from seisflows.config import config, loadpy, tilde_expand, Dict

# run test
if __name__ == '__main__':

    parameters = Dict(loadpy('parameters.py'))
    paths = Dict(tilde_expand(loadpy('paths.py')))

    # register parameters and paths
    sys.modules['seisflows_parameters'] = parameters
    sys.modules['seisflows_paths'] = paths

    # instantiate and register objects
    config()

    # create handles
    workflow = sys.modules['seisflows_workflow']
    system = sys.modules['seisflows_system']

    #cd(args.workdir)

    # submit workflow
    system.submit(workflow)

